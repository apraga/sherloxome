//! # Insilico Controls
//! Generate control variants by sampling clinvar data
use crate::check_deps;
use crate::dbsnp;
use crate::download_blocking;
use crate::run::{Run, run_from_filename};
use crate::setup::{SamplesheetRow, silico_row};
use crate::simuscop;
use crate::varben;
use log;
use noodles::bed;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::io::Write as VCFWrite;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::info::field::Value;
use noodles::vcf::variant::record_buf::AlternateBases as AltBasesBuf;
use rand::RngExt;
use rand::prelude::IteratorRandom;
use serde::Deserialize;
use std::collections::{BTreeSet, HashMap};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::thread;

/// [silico.simuscop] — presence enables simuscop FASTQ generation
/// dbSNP data is always enabled, see ([silico.simuscop.dbsnp]) and the [dbSNP setup
/// guide](https://apraga.github.io/sherloxome/022-dbsnp.html#configuration).
#[derive(Deserialize, Debug)]
pub struct SilicoSimuscopConfig {
    pub capture: String,
    /// Path to a pre-built seqToProfile profile directory. See the [GIAB
    /// example](https://apraga.github.io/sherloxome/0211-simuscop.html#giab-example) for how to
    /// obtain one.
    pub profile: PathBuf,
    /// Target mean sequencing coverage over the capture region. simuscop's own `coverage`
    /// parameter behaves more like a peak/max than a realized mean (see
    /// `MEAN_COVERAGE_REALIZATION`), so this value is scaled up before being written to the
    /// simuReads config.
    pub coverage: u32,
}

/// [silico.varben] — presence enables varben BAM editing
#[derive(Deserialize, Debug)]
pub struct SilicoVarbenConfig {
    pub capture: String,
    /// BAM files to edit, one per patient/sequencer/depth combination to cover. Each must
    /// already exist locally (no URL support — the filename itself is how the run's
    /// patient/sequencer/depth are recovered, via the SAMPLE_SEQUENCER_CAPTURE_DEPTHx
    /// filenaming scheme) and parse to the same capture kit as [silico] `capture`.
    pub bam_files: Vec<PathBuf>,
    /// Minimum read depth required to edit a position (--mindepth)
    pub mindepth: Option<u32>,
}

#[derive(Deserialize, Debug)]
pub struct SilicoConfig {
    pub capture: String,
    /// Local ClinVar VCF path; if absent the file is downloaded from NCBI.
    pub clinvar: Option<PathBuf>,
    /// Number of clinvar variants to sample to insert in the BAM
    pub nb_variants: Option<u32>,
    /// Output directory for intermediate files; defaults to "data/exp_raw".
    pub outdir: Option<PathBuf>,
    /// Simuscop FASTQ generation config ([silico.simuscop]); absence disables it.
    pub simuscop: Option<SilicoSimuscopConfig>,
    /// Varben BAM editing config ([silico.varben]); absence disables it.
    pub varben: Option<SilicoVarbenConfig>,
}

impl SilicoSimuscopConfig {
    /// Build the config from a profile filename alone: capture and coverage are the ones
    /// encoded in `SEQUENCER_CAPTURE_DEPTHx.profile`.
    pub fn from_profile(profile: &Path) -> Result<Self, Box<dyn Error>> {
        let run = simuscop::run_from_profile(profile)?;
        Ok(Self {
            capture: run.capture,
            profile: profile.to_path_buf(),
            coverage: run.depth,
        })
    }
}

/// Where the silico intermediate files and FASTQ are written
pub fn silico_outdir(silico: &SilicoConfig) -> PathBuf {
    silico
        .outdir
        .clone()
        .unwrap_or_else(|| PathBuf::from("data/exp_raw"))
}

/// Sampled clinvar variants, shared by all the runs (see [`prepare_controls`])
fn sampled_clinvar_path(outdir: &Path, clinvar_capture: &str) -> PathBuf {
    outdir.join(format!("clinvar_{clinvar_capture}.vcf.gz"))
}

/// dbSNP variants given to simuReads for a capture kit (see [`prepare_controls`])
fn simuscop_snp_path(outdir: &Path, capture: &str) -> PathBuf {
    outdir.join(format!("dbsnp_{capture}.snp"))
}

/// Write the inputs shared by all the runs, must be done before generating FASTQ
/// - sample clinvar variants inside `bed` (or reuse the ones already sampled)
/// - for every capture of `simuscop_captures`, sample dbSNP variants
/// - if there are varben `bams`, make sure the BWA index of `fasta` exists
pub fn prepare_controls(
    silico: &SilicoConfig,
    bed: PathBuf,
    clinvar_capture: &str,
    fasta: &Path,
    simuscop_captures: &[String],
    bams: &[PathBuf],
) -> Result<(), Box<dyn Error>> {
    let mut tools = vec!["tabix"];
    if !simuscop_captures.is_empty() {
        tools.push("bcftools");
    }
    check_deps(tools);
    let outdir = silico_outdir(silico);
    let clinvar_vcf = sampled_clinvar_path(&outdir, clinvar_capture);
    sample_clinvar(
        silico.clinvar.clone(),
        bed,
        50,
        silico.nb_variants,
        clinvar_vcf.clone(),
    )?;

    for capture in simuscop_captures.iter().collect::<BTreeSet<_>>() {
        let dbsnp_vcf = dbsnp::sample_dbsnp(capture, &clinvar_vcf, &outdir)?;
        dbsnp::write_snp_input(&dbsnp_vcf, &simuscop_snp_path(&outdir, capture))?;
    }
    if !bams.is_empty() {
        varben::ensure_bwa_index(&fasta.to_path_buf())?;
    }
    Ok(())
}

/// Variants sampled by [`prepare_controls`]
fn load_sampled_clinvar(
    outdir: &Path,
    clinvar_capture: &str,
) -> Result<(Vec<RecordBuf>, vcf::Header), Box<dyn Error>> {
    let path = sampled_clinvar_path(outdir, clinvar_capture);
    if !path.exists() {
        return Err(format!(
            "{} not found: run `sherloxome setup prepare` first",
            path.display()
        )
        .into());
    }
    read_vcf(&path)
}

/// Check that a varben BAM exists and follows the filenaming scheme
fn varben_run(bam: &PathBuf) -> Result<Run, Box<dyn Error>> {
    if !bam.exists() {
        return Err(format!("BAM file not found: {}", bam.display()).into());
    }
    run_from_filename(bam).ok_or_else(|| {
        format!(
            "BAM file {:?} does not follow the SAMPLE_SEQUENCER_CAPTURE_DEPTHx filenaming scheme",
            bam
        )
        .into()
    })
}

/// Generate controls from clinvar data and either a BAM file (real patient) or 100% in silico
/// Runs [`prepare_controls`], then every varben and simuscop run of the configuration.
/// Returns a list of samplesheet rows for writing
pub fn generate_controls(
    silico: &SilicoConfig,
    bed: PathBuf,
    clinvar_capture: &str,
    fasta: PathBuf,
) -> Result<Vec<SamplesheetRow>, Box<dyn Error>> {
    let mut tools = vec!["tabix", "bcftools"];
    if silico.simuscop.is_some() {
        tools.push("simuReads");
    }
    if silico.varben.is_some() {
        tools.extend(["bwa", "samtools", "muteditor"]);
    }
    check_deps(tools);

    let captures: Vec<String> = silico.simuscop.iter().map(|s| s.capture.clone()).collect();
    let bams: Vec<PathBuf> = silico
        .varben
        .iter()
        .flat_map(|v| v.bam_files.clone())
        .collect();
    prepare_controls(
        silico,
        bed.clone(),
        clinvar_capture,
        &fasta,
        &captures,
        &bams,
    )?;

    let mut rows: Vec<SamplesheetRow> = Vec::new();
    if let Some(varben) = &silico.varben {
        try_generate_varben(silico, clinvar_capture, &fasta, varben, &mut rows)?
    }
    if let Some(simuscop) = &silico.simuscop {
        rows.push(generate_simuscop(
            silico,
            clinvar_capture,
            &bed,
            &fasta,
            simuscop,
        )?);
    }
    Ok(rows)
}

fn try_generate_varben(
    silico: &SilicoConfig,
    clinvar_capture: &str,
    fasta: &PathBuf,
    varben: &SilicoVarbenConfig,
    rows: &mut Vec<SamplesheetRow>,
) -> Result<(), Box<dyn Error>> {
    if varben.bam_files.is_empty() {
        log::error!("[silico.varben] requires at least one entry in bam_files");
        return Err("[silico.varben] requires at least one entry in bam_files".into());
    }
    for bam in &varben.bam_files {
        let run = varben_run(bam)?;
        if run.capture != varben.capture {
            return Err(format!(
                    "BAM {:?} is for capture '{}' but [silico.varben] capture is '{}': list only BAMs for \
                     this capture kit here, and run `setup` once per kit to cover several",
                    bam, run.capture, varben.capture
                )
                .into());
        }
        rows.push(generate_varben_single(silico, clinvar_capture, fasta, bam)?);
    }
    Ok(())
}

/// Run varben to insert sampled clinvar variants into a single `bam` and convert it to FASTQ.
/// [`prepare_controls`] must have been run before.
pub fn generate_varben_single(
    silico: &SilicoConfig,
    clinvar_capture: &str,
    fasta: &PathBuf,
    bam: &PathBuf,
) -> Result<SamplesheetRow, Box<dyn Error>> {
    check_deps(vec!["bwa", "samtools", "tabix", "muteditor"]);
    let run = varben_run(bam)?;
    let outdir = silico_outdir(silico);
    let mindepth = silico.varben.as_ref().and_then(|v| v.mindepth);
    let (variants, header) = load_sampled_clinvar(&outdir, clinvar_capture)?;

    let (fq1, fq2) = varben::generate_controls_bam(
        bam,
        &run.capture,
        fasta,
        &variants,
        &header,
        mindepth,
        &outdir,
    )?;
    Ok(silico_row("varben", &run.capture, fq1, fq2))
}

/// One simuscop run: generate a FASTQ from `simuscop.profile`. Needs [`prepare_controls`] to
/// have run for `simuscop.capture`.
pub fn generate_simuscop(
    silico: &SilicoConfig,
    clinvar_capture: &str,
    bed: &PathBuf,
    fasta: &PathBuf,
    simuscop: &SilicoSimuscopConfig,
) -> Result<SamplesheetRow, Box<dyn Error>> {
    check_deps(vec!["simuReads", "tabix"]);
    // simuscop's `coverage` parameter behaves like a peak/max rather than a realized mean:
    // This factor was measured empirically
    const MEAN_COVERAGE_REALIZATION: f64 = 0.65;
    let coverage = (simuscop.coverage as f64 / MEAN_COVERAGE_REALIZATION).round() as u32;
    log::debug!(
        "Scaling coverage {} -> {coverage} to convert to mean",
        simuscop.coverage
    );

    let outdir = silico_outdir(silico);
    let (variants, header) = load_sampled_clinvar(&outdir, clinvar_capture)?;
    let capture = simuscop.capture.as_str();
    let snp_path = simuscop_snp_path(&outdir, capture);
    if !snp_path.exists() {
        return Err(format!(
            "{} not found: run `sherloxome setup prepare` first",
            snp_path.display()
        )
        .into());
    }

    let (fq1, fq2) = simuscop::generate_controls_fastq(
        bed,
        fasta,
        &simuscop.profile,
        &variants,
        &header,
        &outdir,
        coverage,
        snp_path,
    )?;
    Ok(silico_row("simuscop", capture, fq1, fq2))
}

/// Rebuild the samplesheet rows of every silico FASTQ pair found in `outdir`: the filename
/// gives the capture kit and the tool (`..._simuscop_1.fq.gz`, `..._varben_1.fq.gz`).
/// Rows are sorted by sample name.
pub fn silico_rows_from_disk(outdir: &Path) -> Result<Vec<SamplesheetRow>, Box<dyn Error>> {
    let mut rows = Vec::new();
    for entry in
        std::fs::read_dir(outdir).map_err(|e| format!("Cannot read {}: {e}", outdir.display()))?
    {
        let fq1 = entry?.path();
        let Some(name) = fq1.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        let Some(prefix) = name.strip_suffix("_1.fq.gz") else {
            continue;
        };
        let Some(run) = run_from_filename(&PathBuf::from(name)) else {
            continue;
        };
        let Some(tool) = run
            .silico
            .as_deref()
            .filter(|s| matches!(*s, "simuscop" | "varben"))
        else {
            continue;
        };
        let fq2 = outdir.join(format!("{prefix}_2.fq.gz"));
        if !fq2.exists() {
            log::warn!("Skipping {:?}: {:?} is missing", fq1, fq2);
            continue;
        }
        rows.push(silico_row(tool, &run.capture, fq1, fq2));
    }
    rows.sort_by(|a, b| a.sample.cmp(&b.sample));
    Ok(rows)
}

/// Generate controls from a BAM file and returns 2 fastq
fn download_clinvar() -> PathBuf {
    log::debug!("Downloading clinvar vcf...");
    let mut url: String =
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz".to_string();
    let vcf = PathBuf::from("data/exp_raw/clinvar.vcf.gz");
    download_blocking(&url, &vcf);

    url.push_str(".tbi");
    let tbi = PathBuf::from("data/exp_raw/clinvar.vcf.gz.tbi");
    download_blocking(&url, &tbi);
    vcf
}

/// Pseudo-interval tree..
type BedIndex = HashMap<String, Vec<(usize, usize)>>;

/// For each chromosome, store sorted list of (start, end) 0-based half-open intervals
/// Assume bed is stored
fn load_bed(bed: &PathBuf) -> Result<BedIndex, Box<dyn Error>> {
    let mut reader = File::open(bed)
        .map(BufReader::new)
        .map(bed::io::Reader::<3, _>::new)?;

    let mut index = BedIndex::new();
    let mut record = bed::Record::default();

    while reader.read_record(&mut record)? != 0 {
        let chrom_raw = std::str::from_utf8(record.reference_sequence_name())?;
        let chrom = chrom_raw
            .strip_prefix("chr")
            .unwrap_or(chrom_raw)
            .to_string();
        let start = usize::from(record.feature_start()?);
        let end = usize::from(record.feature_end().ok_or("missing BED end position")??);
        index.entry(chrom).or_default().push((start, end));
    }
    Ok(index)
}

/// Is a 1-based VCF position inside any BED interval on this chromosome?
fn in_capture(bed: &BedIndex, chrom: &str, pos: u64) -> bool {
    // Convert 1-based VCF pos to 0-based
    let pos0 = pos as usize - 1;

    let Some(intervals) = bed.get(chrom) else {
        return false;
    };

    // Binary search: find the last interval that starts at or before pos0
    let i = intervals.partition_point(|(start, _)| *start <= pos0);

    // Check if pos0 falls before the end of that interval
    i.checked_sub(1)
        .map(|j| intervals[j].1 > pos0)
        .unwrap_or(false)
}

/// Pathogenic variant are preferred
fn clnsig_priority(sig: &str) -> Option<u8> {
    match sig {
        "Pathogenic" => Some(0),
        "Likely_pathogenic" => Some(1),
        "Uncertain_significance" => Some(2),
        _ => None,
    }
}
/// Only keep VOUS, Pathogenic or Likely Pathogenic or VOUS variant (whatever the number of submisson)
/// If there are several variants, only keep the most severe
/// Info field is complicated... Some(Ok(Some(Array([Ok(Some("Uncertain_significance"))]))))
fn is_not_benign(info: &vcf::record::Info, header: &vcf::Header) -> bool {
    info.get(header, "CLNSIG")
        .and_then(|r| r.ok())
        .flatten()
        .is_some_and(|v| match v {
            Value::String(s) => clnsig_priority(s.as_ref()).is_some(),
            Value::Array(arr) => match arr {
                noodles::vcf::variant::record::info::field::value::Array::String(arr) => arr
                    .iter()
                    .filter_map(|s| s.ok().flatten())
                    .any(|s| clnsig_priority(s.as_ref()).is_some()),
                _ => false,
            },
            _ => false,
        })
}

/// Only keep SNV
fn is_snv(record: &vcf::Record) -> bool {
    record.reference_bases().len() == 1
        && record
            .alternate_bases()
            .iter()
            .filter_map(|a| a.ok())
            .all(|a| a.len() == 1)
}
/// Sort chromsome by natural ornder
fn chrom_order(chrom: &str) -> (u8, u32) {
    match chrom.trim_start_matches("chr") {
        "X" => (1, 0),
        "Y" => (2, 0),
        "M" | "MT" => (3, 0),
        n => (0, n.parse().unwrap_or(u32::MAX)),
    }
}

/// Sample clinvar variant randomly.
/// Prefix "chr" to chromosome names
fn sample_clinvar_variants(
    reader: &mut vcf::io::Reader<bgzf::io::Reader<File>>,
    header: &vcf::Header,
    capture: &BedIndex,
    spacing: u32,
    n: u32,
    rng: &mut impl RngExt,
) -> Vec<RecordBuf> {
    let mut last_pos: HashMap<String, u64> = HashMap::new();
    let mut selected: Vec<RecordBuf> = reader
        .records()
        .filter_map(|r| r.ok())
        // Filter first as we need INFO field
        .filter(|record| keep_variant(record, header, capture, &mut last_pos, spacing))
        // The add chr prefilx
        .map(|record| add_chr_prefix(&record))
        .sample(
            rng,
            n.try_into().expect("Fails to convert nb variants to usize"),
        );
    sort_by_chromosome(&mut selected);
    log::debug!("Selected {} variants", selected.len());
    selected
}

pub fn sort_by_chromosome(variants: &mut Vec<RecordBuf>) {
    variants.sort_by(|a, b| {
        chrom_order(a.reference_sequence_name())
            .cmp(&chrom_order(b.reference_sequence_name()))
            .then_with(|| {
                let ap = a.variant_start().map(usize::from).unwrap_or(0);
                let bp = b.variant_start().map(usize::from).unwrap_or(0);
                ap.cmp(&bp)
            })
    });
}

/// Tabix-index a bgzipped VCF (writes `{vcf}.tbi`), overwriting any existing index.
pub fn index_vcf(vcf: &Path) -> Result<(), Box<dyn Error>> {
    log::debug!("Indexing {:?}", vcf);
    let status = Command::new("tabix")
        .args(["-f", "-p", "vcf"])
        .arg(vcf)
        .status()?;
    if !status.success() {
        return Err(format!("tabix exited with status {status}").into());
    }
    Ok(())
}

/// Select n clinvar pathogenic SNV inside the capture kit 50bp apart
/// Output is a list of variant. Variants are also written in a data/exp_raw/clinvar_$CAPTURE.vcf
/// The most efficient way is to parse clinvar once to get variant in the capture kit and 50bp apart.
/// In a second pass, sample randomly n of them.
pub fn sample_clinvar(
    clinvar_vcf: Option<PathBuf>,
    bed: PathBuf,
    spacing: u32,
    nb_variants: Option<u32>,
    vcf_out: PathBuf,
) -> Result<(Vec<RecordBuf>, vcf::Header), Box<dyn Error>> {
    let clinvar_path = match clinvar_vcf {
        Some(p) => p,
        None => download_clinvar(),
    };

    let n = match nb_variants {
        Some(n) => n,
        None => 1000,
    };

    let capture = load_bed(&bed)?;
    let mut reader = File::open(&clinvar_path)
        .map(bgzf::io::Reader::new)
        .map(vcf::io::Reader::new)
        .expect("Failed to open clinvar file");
    let header = reader.read_header()?;

    if vcf_out.exists() {
        log::debug!(
            "Skip sampling clinvar as output files already exists: {:?}",
            vcf_out
        );
        let (variants, header) = read_vcf(&vcf_out)?;
        Ok((variants, header))
    } else {
        log::debug!("Sampling {n} variants for insertion");
        let mut rng = rand::rng();
        let variants =
            sample_clinvar_variants(&mut reader, &header, &capture, spacing, n, &mut rng);

        write_sampled_clinvar_vcf(&variants, &header, &vcf_out)?;
        log::debug!("Wrote {:?}", vcf_out);
        Ok((variants, header))
    }
}

/// Read a bgzf-compressed VCF and return all records and the header
fn read_vcf(path: &PathBuf) -> Result<(Vec<RecordBuf>, vcf::Header), Box<dyn Error>> {
    let mut reader = File::open(path)
        .map(bgzf::io::Reader::new)
        .map(vcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let records = reader.record_bufs(&header).filter_map(|r| r.ok()).collect();
    Ok((records, header))
}

/// Write sampled clinvar as VCF and index it.
/// Assume there is no chr prefix for chromosome
fn write_sampled_clinvar_vcf(
    variants: &Vec<RecordBuf>,
    header: &vcf::Header,
    vcf_out: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let mut writer = File::create(&vcf_out)
        .map(bgzf::io::Writer::new)
        .map(vcf::io::Writer::new)?;
    writer.write_header(&header)?;

    for record in variants {
        writer.write_variant_record(header, record)?;
    }
    // Drop first: the inner bgzf writer only flushes its last block and writes the BGZF EOF
    // marker on Drop, and tabix needs that on disk before it can index the file.
    drop(writer);
    index_vcf(vcf_out)
}

/// We cannot edit a record directly, recreate it without INFO and sample.
/// This is just for clinvar data for quick check
fn add_chr_prefix(record: &vcf::Record) -> RecordBuf {
    let chrom_raw = record.reference_sequence_name();
    let chrom_nb = chrom_raw.strip_prefix("chr").unwrap_or(chrom_raw);

    let alts: Vec<String> = record
        .alternate_bases()
        .iter()
        .filter_map(|a| a.ok())
        .map(|a| a.to_string())
        .collect();

    let mut builder = RecordBuf::builder()
        .set_reference_sequence_name(format!("chr{}", chrom_nb))
        .set_reference_bases(record.reference_bases().to_string())
        .set_alternate_bases(AltBasesBuf::from(alts));

    if let Some(Ok(pos)) = record.variant_start() {
        builder = builder.set_variant_start(pos);
    }

    builder.build()
}

/// For varian in the capture file and SNV and not bening, save it for writing
/// Update lats saved position (`last_pos)`
/// Input VCF may use chr prefix
fn keep_variant(
    record: &vcf::Record,
    header: &vcf::Header,
    capture: &BedIndex,
    last_pos: &mut HashMap<String, u64>,
    spacing: u32,
) -> bool {
    let Some(Ok(start)) = record.variant_start() else {
        return false;
    };
    let pos = usize::from(start) as u64;
    let chrom_raw = record.reference_sequence_name();
    let chrom = chrom_raw.strip_prefix("chr").unwrap_or(chrom_raw);
    let last = last_pos.get(chrom).copied().unwrap_or(0);

    let ok = in_capture(capture, chrom, pos)
        && is_not_benign(&record.info(), header)
        && is_snv(record)
        && pos != last
        && pos.abs_diff(last) >= spacing.into();

    if ok {
        last_pos.insert(chrom.to_string(), pos);
    }
    ok
}

pub fn nb_threads() -> usize {
    thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn simuscop_config_has_required_fields() {
        let dir = std::env::temp_dir().join("simuscop_config_test");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("test.conf");

        simuscop::write_config(
            &config_path,
            Path::new("/ref/genome.fa"),
            Path::new("/data/sample.profile"),
            Path::new("/data/clinvar.simuscop"),
            Path::new("/data/capture.bed"),
            "agilent-col6a1",
            Path::new("/data/simuscop_out"),
            50,
            PathBuf::from("snp.txt"),
        )
        .unwrap();

        let content = fs::read_to_string(&config_path).unwrap();
        assert!(content.contains("ref = /ref/genome.fa"), "missing ref");
        assert!(
            content.contains("profile = /data/sample.profile"),
            "missing profile"
        );
        assert!(
            content.contains("variation = /data/clinvar.simuscop"),
            "missing variation"
        );
        assert!(
            content.contains("target = /data/capture.bed"),
            "missing target"
        );
        assert!(content.contains("name = agilent-col6a1"), "missing name");
        assert!(
            content.contains("output = /data/simuscop_out"),
            "missing output"
        );
        assert!(content.contains("layout = PE"), "missing layout");
        assert!(content.contains("coverage = 50"), "missing coverage");
        assert!(content.contains("threads = "), "missing threads");
    }

    #[test]
    fn simuscop_config_coverage_matches_input() {
        let dir = std::env::temp_dir().join("simuscop_coverage_test");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("test.conf");

        simuscop::write_config(
            &config_path,
            Path::new("/ref/genome.fa"),
            Path::new("/profile"),
            Path::new("/variation"),
            Path::new("/bed"),
            "sample",
            Path::new("/out"),
            100,
            PathBuf::from("snp.txt"),
        )
        .unwrap();

        let content = fs::read_to_string(&config_path).unwrap();
        assert!(content.contains("coverage = 100"));
    }

    #[test]
    fn simuscop_config_includes_snp_when_set() {
        let dir = std::env::temp_dir().join("simuscop_snp_test");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("test.conf");

        simuscop::write_config(
            &config_path,
            Path::new("/ref/genome.fa"),
            Path::new("/profile"),
            Path::new("/variation"),
            Path::new("/bed"),
            "sample",
            Path::new("/out"),
            50,
            PathBuf::from("/data/dbsnp_agilent-col6a1.snp"),
        )
        .unwrap();

        let content = fs::read_to_string(&config_path).unwrap();
        assert!(
            content.contains("snp = /data/dbsnp_agilent-col6a1.snp"),
            "missing snp"
        );
    }

    #[test]
    fn simuscop_config_from_profile_filename() {
        let conf = SilicoSimuscopConfig::from_profile(Path::new(
            "data/ref/profiles/novaseq_idt_100x.profile",
        ))
        .unwrap();
        assert_eq!(conf.capture, "idt");
        assert_eq!(conf.coverage, 100);
        assert_eq!(
            conf.profile,
            PathBuf::from("data/ref/profiles/novaseq_idt_100x.profile")
        );
    }

    #[test]
    fn simuscop_config_from_profile_keeps_dash_in_capture() {
        let conf =
            SilicoSimuscopConfig::from_profile(Path::new("hiseq4000_agilent-col6a1_50x.profile"))
                .unwrap();
        assert_eq!(conf.capture, "agilent-col6a1");
        assert_eq!(conf.coverage, 50);
    }

    #[test]
    fn simuscop_config_from_bad_profile_name_fails() {
        assert!(SilicoSimuscopConfig::from_profile(Path::new("mine.profile")).is_err());
    }

    #[test]
    fn runs_refuse_to_sample_when_not_prepared() {
        let dir = std::env::temp_dir().join("silico_not_prepared_test");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();

        let err = load_sampled_clinvar(&dir, "agilent-idt-truseq")
            .unwrap_err()
            .to_string();
        assert!(err.contains("setup prepare"), "{err}");
        // Nothing was sampled in place of the missing file
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 0);
    }

    #[test]
    fn silico_rows_from_disk_keeps_complete_silico_pairs() {
        let dir = std::env::temp_dir().join("silico_rows_from_disk_test");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        for f in [
            "nopatient_novaseq_idt_50x_simuscop_1.fq.gz",
            "nopatient_novaseq_idt_50x_simuscop_2.fq.gz",
            "HG002_hiseq4000_agilent_50x_varben_1.fq.gz",
            "HG002_hiseq4000_agilent_50x_varben_2.fq.gz",
            // no mate
            "nopatient_novaseq_truseq_50x_simuscop_1.fq.gz",
            // real patient, not silico
            "HG002_novaseq_idt_50x_1.fq.gz",
            "HG002_novaseq_idt_50x_2.fq.gz",
            // not a fastq
            "nopatient_novaseq_idt_50x_simuscop.conf",
        ] {
            fs::write(dir.join(f), "").unwrap();
        }

        let rows = silico_rows_from_disk(&dir).unwrap();
        let got: Vec<(&str, &str, &str)> = rows
            .iter()
            .map(|r| (r.patient.as_str(), r.capture.as_str(), r.sample.as_str()))
            .collect();
        assert_eq!(
            got,
            vec![
                (
                    "silico-varben",
                    "agilent",
                    "HG002_hiseq4000_agilent_50x_varben"
                ),
                (
                    "silico-simuscop",
                    "idt",
                    "nopatient_novaseq_idt_50x_simuscop"
                ),
            ]
        );
        assert!(
            rows[1]
                .fastq_2
                .ends_with("nopatient_novaseq_idt_50x_simuscop_2.fq.gz")
        );
    }
}
