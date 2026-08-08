//! # Insilico Controls
//! Generate control variants by sampling clinvar data
use crate::download_blocking;
use crate::setup::{SamplesheetRow, silico_row};
use crate::simuscop;
use crate::varben;
use crate::{check_deps, resolve_bam};
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
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::thread;

/// [silico.simuscop] — presence enables simuscop FASTQ generation
#[derive(Deserialize, Debug)]
pub struct SilicoSimuscopConfig {
    /// Path to a pre-built seqToProfile profile directory. Mutually exclusive with `vcf`.
    pub profile: Option<PathBuf>,
    /// VCF of germline variants called from bam_file (e.g. via GATK HaplotypeCaller).
    /// Required when `profile` is absent; seqToProfile is run to build the profile.
    pub vcf: Option<PathBuf>,
    /// Sequencing coverage
    pub coverage: u32,
}

/// [silico.varben] — presence enables varben BAM editing
#[derive(Deserialize, Debug)]
pub struct SilicoVarbenConfig {
    /// Minimum read depth required to edit a position (--mindepth)
    pub mindepth: Option<u32>,
}

#[derive(Deserialize, Debug)]
pub struct SilicoConfig {
    pub capture: String,
    /// A BAM file is required for varben but not nor simuscop
    pub bam_file: Option<String>,
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

/// Generate controls from clinvar data and either a BAM file (real patient) or 100% in silico
/// Returns a list of samplesheet rows for writing
pub fn generate_controls(
    silico: &SilicoConfig,
    bed: PathBuf,
    capture: &str,
    fasta: PathBuf,
) -> Result<Vec<SamplesheetRow>, Box<dyn Error>> {
    check_deps(&["bwa", "samtools", "muteditor", "seqToProfile"]);

    let outdir = silico
        .outdir
        .clone()
        .unwrap_or_else(|| PathBuf::from("data/exp_raw"));

    let mut rows: Vec<SamplesheetRow> = Vec::new();

    let vcf_out = outdir.join(format!("clinvar_{capture}.vcf.gz"));
    let (variants, header) = sample_clinvar(
        silico.clinvar.clone(),
        PathBuf::from(bed.clone()),
        50,
        silico.nb_variants,
        vcf_out,
    )?;
    if let Some(varben) = &silico.varben {
        generate_controls_varben(
            &silico, &bed, capture, &fasta, &variants, &header, &outdir, varben, &mut rows,
        )?;
    }
    if let Some(simuscop) = &silico.simuscop {
        generate_controls_simuscop(
            &silico, &bed, capture, &fasta, &variants, &header, &outdir, simuscop, &mut rows,
        )?;
    }
    Ok(rows)
}

fn generate_controls_varben(
    silico: &SilicoConfig,
    bed: &PathBuf,
    capture: &str,
    fasta: &PathBuf,
    variants: &Vec<RecordBuf>,
    header: &vcf::Header,
    outdir: &PathBuf,
    varben: &SilicoVarbenConfig,
    rows: &mut Vec<SamplesheetRow>,
) -> Result<(), Box<dyn Error>> {
    if let Some(bam_file) = &silico.bam_file {
        let bam_path = resolve_bam(bam_file, &outdir)?;
        let (fq1, fq2) = varben::generate_controls_bam(
            &bam_path,
            bed,
            capture,
            fasta,
            variants,
            header,
            varben.mindepth,
            outdir,
        )?;
        rows.push(silico_row("varben", fq1, fq2));
        Ok(())
    } else {
        log::error!("[silico.varben] requires a BAM file");
        return Err("[silico.varben] requires a BAM file".into());
    }
}

fn generate_controls_simuscop(
    silico: &SilicoConfig,
    bed: &PathBuf,
    capture: &str,
    fasta: &PathBuf,
    variants: &Vec<RecordBuf>,
    header: &vcf::Header,
    outdir: &PathBuf,
    simuscop: &SilicoSimuscopConfig,
    rows: &mut Vec<SamplesheetRow>,
) -> Result<(), Box<dyn Error>> {
    let (fq1, fq2) = simuscop::generate_controls_fastq(
        &silico.bam_file,
        &bed,
        &capture,
        &fasta,
        simuscop.vcf.as_ref(),
        simuscop.profile.as_ref(),
        &variants,
        &outdir,
        simuscop.coverage,
    )?;
    rows.push(silico_row("simuscop", fq1, fq2));
    Ok(())
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
    Ok(())
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
        )
        .unwrap();

        let content = fs::read_to_string(&config_path).unwrap();
        assert!(content.contains("coverage = 100"));
    }
}
