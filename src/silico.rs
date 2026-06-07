//! # Insilico Controls
//! Generate control variants by sampling clinvar data
use crate::setup::{SamplesheetRow, silico_row};
use crate::{check_deps, resolve_bam};
use crate::{download_blocking, ref_dir};
use log;
use noodles::bed;
use noodles::bgzf;
use noodles::fasta;
use noodles::vcf;
use noodles::vcf::header::record::value::{Map, map::Format};
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::io::Write as VCFWrite;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::info::field::Value;
use noodles::vcf::variant::record::samples::keys::key as sample_key;
use noodles::vcf::variant::record_buf::AlternateBases as AltBasesBuf;
use noodles::vcf::variant::record_buf::samples::{
    Keys as SampleKeys, Samples as SamplesBuf, sample::Value as SampleValue,
};
use rand::RngExt;
use rand::prelude::IteratorRandom;
use serde::Deserialize;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitStatus;
use std::process::{Command, Stdio};
use std::thread;

/// [silico.simuscop] — presence enables simuscop FASTQ generation
#[derive(Deserialize, Debug)]
pub struct SilicoSimuscopConfig {
    /// VCF of germline variants called from bam_file (e.g. via GATK HaplotypeCaller)
    pub vcf: PathBuf,
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
    /// A BAM file is always required. Varben will modify it, simuscop will define a sequencing profile from it.
    pub bam_file: String,
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
    let bam_path = resolve_bam(&silico.bam_file, &outdir)?;

    if let Some(varben) = &silico.varben {
        let (fq1, fq2) = generate_controls_bam(
            &bam_path, &bed, &capture, &fasta, &variants, header, varben.mindepth, &outdir,
        )?;
        rows.push(silico_row("varben", fq1, fq2));
    }
    if let Some(simuscop) = &silico.simuscop {
        let (fq1, fq2) = generate_controls_fastq(
            &bam_path, &bed, &capture, &fasta, &simuscop.vcf, &variants, &outdir, simuscop.coverage,
        )?;
        rows.push(silico_row("simuscop", fq1, fq2));
    }
    Ok(rows)
}

/// Generate controls into a in-silico FASTQ. The BAM file is needed to generate a sequencing profile
fn generate_controls_fastq(
    bam: &PathBuf,
    bed: &PathBuf,
    capture: &str,
    fasta: &PathBuf,
    vcf: &PathBuf,
    variants: &Vec<RecordBuf>,
    outdir: &PathBuf,
    coverage: u32,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let profile_dir = ensure_profile(bam, vcf, fasta, bed, outdir)?;

    let variation = outdir.join(format!("clinvar_{capture}.simuscop"));
    write_simuscop_input(variants, &variation, capture)?;

    let simuscop_outdir = outdir.join(format!("simuscop_{capture}"));
    std::fs::create_dir_all(&simuscop_outdir)?;

    let config_path = outdir.join(format!("simuscop_{capture}.conf"));
    write_simuscop_config(&config_path, fasta, &profile_dir, &variation, bed, capture, &simuscop_outdir, coverage)?;

    run_simu_reads(&config_path, capture, &simuscop_outdir)
}

/// Write a simuReads config file
pub fn write_simuscop_config(
    config_path: &Path,
    fasta: &Path,
    profile_dir: &Path,
    variation: &Path,
    bed: &Path,
    name: &str,
    output_dir: &Path,
    coverage: u32,
) -> Result<(), Box<dyn Error>> {
    let threads = nb_threads();
    let mut w = BufWriter::new(File::create(config_path)?);
    writeln!(w, "ref = {}", fasta.display())?;
    writeln!(w, "profile = {}", profile_dir.display())?;
    writeln!(w, "variation = {}", variation.display())?;
    writeln!(w, "target = {}", bed.display())?;
    writeln!(w, "name = {}", name)?;
    writeln!(w, "output = {}", output_dir.display())?;
    writeln!(w, "layout = PE")?;
    writeln!(w, "threads = {}", threads)?;
    writeln!(w, "coverage = {}", coverage)?;
    w.flush()?;
    Ok(())
}

/// Run simuReads with the given config. Returns (fq1, fq2) output paths.
fn run_simu_reads(
    config_path: &Path,
    name: &str,
    output_dir: &Path,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let status = Command::new("simuReads")
        .arg(config_path.to_str().ok_or("Invalid simuReads config path")?)
        .status()?;

    if !status.success() {
        return Err(format!("simuReads exited with status {status}").into());
    }

    let fq1 = output_dir.join(format!("{name}_1.fq.gz"));
    let fq2 = output_dir.join(format!("{name}_2.fq.gz"));
    Ok((fq1, fq2))
}

/// Build a seqToProfile sequencing profile from a normal BAM.
/// The profile directory is derived from the BAM stem and reused on subsequent runs.
fn ensure_profile(
    bam: &PathBuf,
    vcf: &PathBuf,
    fasta: &PathBuf,
    bed: &PathBuf,
    outdir: &PathBuf,
) -> Result<PathBuf, Box<dyn Error>> {
    let bam_stem = bam.file_stem().unwrap().to_string_lossy();
    let profile_dir = outdir.join(format!("{bam_stem}.profile"));

    if profile_dir.exists() {
        log::debug!("Reusing existing profile: {:?}", profile_dir);
        return Ok(profile_dir);
    }

    std::fs::create_dir_all(&profile_dir)?;
    log::info!("Generating sequencing profile with seqToProfile");

    let status = Command::new("seqToProfile")
        .args([
            "-b",
            bam.to_str().ok_or("Invalid BAM path")?,
            "-v",
            vcf.to_str().ok_or("Invalid VCF path")?,
            "-r",
            fasta.to_str().ok_or("Invalid FASTA path")?,
            "-t",
            bed.to_str().ok_or("Invalid BED path")?,
            "-o",
            profile_dir.to_str().ok_or("Invalid profile dir path")?,
        ])
        .status()?;

    if !status.success() {
        // Remove the empty directory so the next run retries cleanly
        let _ = std::fs::remove_dir(&profile_dir);
        return Err(format!("seqToProfile exited with status {status}").into());
    }

    Ok(profile_dir)
}

/// Generate controls from a BAM file and returns 2 fastq
fn generate_controls_bam(
    bam: &PathBuf,
    bed: &PathBuf,
    capture: &str,
    fasta: &PathBuf,
    variants: &Vec<RecordBuf>,
    header: vcf::Header,
    mindepth: Option<u32>,
    outdir: &PathBuf,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    std::fs::create_dir_all(&outdir)?;

    let bam_stem = bam.file_stem().unwrap().to_string_lossy();
    let bam_parent = bam.parent().unwrap_or(Path::new("."));
    let fq1 = bam_parent.join(format!("{bam_stem}_1.fq.gz"));
    let fq2 = bam_parent.join(format!("{bam_stem}_2.fq.gz"));

    if fq1.exists() && fq2.exists() {
        log::info!(
            "Skipping BAM editing as output fastq already exists {:?}, {:?}",
            fq1,
            fq2
        );
        Ok((fq1, fq2))
    } else {
        log::info!("Editing BAM to insert control");
        let new_bam = edit_bam(&bam, &bed, capture, fasta, variants, header, mindepth)?;
        bam_to_fastq(new_bam, fq1, fq2)
    }
}

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

/// Write VCF record to varben format. Assume chromosoes do NOT have a chr prefix
/// Output generate heterozygous variant (AF betwen 0.4 and 0.6)
fn write_varben(w: &mut impl Write, record: &RecordBuf) -> Result<(), Box<dyn Error>> {
    let mut rng = rand::rng();
    let af: f64 = rng.random_range(0.4..0.6);
    let chrom_raw = record.reference_sequence_name();
    let chrom_nb = chrom_raw.strip_prefix("chr").unwrap_or(chrom_raw);
    let chrom_chr = format!("chr{}", chrom_nb);
    let pos = record.variant_start().unwrap();
    let alt = record.alternate_bases();
    let alt_merged = alt.iter().filter_map(|a| a.ok()).next().unwrap_or(".");
    writeln!(
        w,
        "{}\t{}\t{}\t{}\tsnv\t{}",
        chrom_chr, pos, pos, af, alt_merged
    )?;
    Ok(())
}

/// Write one SNV in simuscop tab-separated variation format (no header line).
/// Always emits `het` since we sample heterozygous AF from clinvar.
fn write_simuscop_variant(
    w: &mut impl Write,
    record: &RecordBuf,
    sample: &str,
) -> Result<(), Box<dyn Error>> {
    let chrom = record.reference_sequence_name();
    let pos = usize::from(record.variant_start().unwrap());
    let ref_base = record.reference_bases().to_lowercase();
    let alt = record
        .alternate_bases()
        .iter()
        .filter_map(|a| a.ok())
        .next()
        .unwrap_or(".");
    writeln!(w, "s\t{sample}\t{chrom}\t{pos}\t{ref_base}\t{alt}\thet")?;
    Ok(())
}

/// Write simuscop input file from a list of variant
fn write_simuscop_input(
    variants: &[RecordBuf],
    out: &Path,
    sample: &str,
) -> Result<(), Box<dyn Error>> {
    let mut w = BufWriter::new(File::create(out)?);
    for v in variants {
        write_simuscop_variant(&mut w, v, sample)?;
    }
    w.flush()?;
    Ok(())
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

fn sort_by_chromosome(variants: &mut Vec<RecordBuf>) {
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

/// Write varben output as VCF successful insertions only, failures are not properly formatted,
pub fn write_varben_as_vcf(
    bam: &PathBuf,
    header: vcf::Header,
    outdir: PathBuf,
    fasta: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let varben_dir = outdir.join("varben");
    let bam_stem = bam.file_stem().unwrap().to_string_lossy();

    let success_list = varben_dir.join("success_list.txt");
    let success_vcf = outdir.join(format!("{bam_stem}_varben.vcf.gz"));
    write_varben_as_vcf_single(success_list, success_vcf, header, fasta)
}

/// Write a mutation file for varben using sampled clinvar data
fn write_varben_input(variants: &Vec<RecordBuf>, mut_out: &PathBuf) -> Result<(), Box<dyn Error>> {
    // Prepare mutation file
    let mut varben_w = BufWriter::new(File::create(&mut_out)?);
    writeln!(varben_w, "#chrom\tstart\tend\tAF\ttype\talt")?;

    for v in variants {
        write_varben(&mut varben_w, v)?;
    }
    varben_w.flush()?;
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

/// Helper function to index a BAM file
pub fn index_bam(bam: &PathBuf) -> Result<(), Box<dyn Error>> {
    let bai = bam.with_extension("bam.bai");
    // Index the BAM
    if !bai.exists() {
        log::debug!("Indexing bam");
        let status = Command::new("samtools")
            .args(["index", &bam.to_string_lossy()])
            .status()
            .expect("Failed to run samtools index");
        if !status.success() {
            return Err(format!("samtools index exited with status {status}").into());
        }
    } else {
        log::debug!("{:?} already exists ", bai);
    }
    Ok(())
}

/// Bam to fastq conversion require hard clips removal
pub fn remove_hard_clips(bam: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    let stem = bam.file_stem().unwrap().to_string_lossy();
    let parent = bam.parent().unwrap_or(Path::new("."));
    let filtered = parent.join(format!("{stem}_nohardclip.bam"));

    if !filtered.exists() {
        log::info!("Removing hard clips");
        let cmd = format!(
            "samtools view -h {} | awk '$6 !~ /H/{{print}}' | samtools view -bS - > {}",
            bam.display(),
            filtered.display()
        );
        let status = Command::new("sh")
            .args(["-c", &cmd])
            .status()
            .expect("Failed to run hard clip filter");

        if !status.success() {
            return Err(format!("samtools exited with status {status}").into());
        }
    } else {
        log::info!("Hardlink removal already done");
    }
    index_bam(&filtered)?;
    Ok(filtered)
}

/// Download and extract the pre-built BWA index from NCBI if absent.
pub fn ensure_bwa_index(fasta: &PathBuf) -> Result<(), Box<dyn Error>> {
    let bwt = PathBuf::from(format!("{}.bwt", fasta.display()));
    if bwt.exists() {
        log::debug!("BWA index already exists");
        return Ok(());
    }
    let tar_name = "GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz";
    let root_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines";
    let url = format!("{root_url}/{tar_name}");
    let tar_path = ref_dir().join(tar_name);
    log::debug!("Downloading BWA index...");
    download_blocking(&url, &tar_path);
    log::debug!("Extracting BWA index...");
    let status = Command::new("tar")
        .args([
            "-xzf",
            tar_path.to_str().ok_or("Invalid tar path")?,
            "-C",
            ref_dir().to_str().ok_or("Invalid ref dir")?,
        ])
        .status()?;
    if !status.success() {
        // Archive is likely corrupt (e.g. interrupted download); remove it so the next run re-downloads.
        let _ = std::fs::remove_file(&tar_path);
        return Err(format!(
            "Failed to extract BWA index archive (tar exited with {status}). \
             The archive may have been corrupt or incompletely downloaded — \
             it has been deleted. Re-run setup to download it again."
        )
        .into());
    }
    Ok(())
}
/// Generate controls from clinvar variant into a BAM file
/// 1. Sample clinvar randomly
/// 2. Generate mutation files from those clinvar variant for varben
/// 3. Apply varben on a single bam file
///
/// If `clinvar_vcf` is None the file is downloaded from NCBI and stored in `outdir`.
pub fn edit_bam(
    bam: &PathBuf,
    bed: &PathBuf,
    capture: &str,
    fasta: &PathBuf,
    variants: &Vec<RecordBuf>,
    header: vcf::Header,
    mindepth: Option<u32>,
) -> Result<PathBuf, Box<dyn Error>> {
    let outdir = PathBuf::from("data/exp_raw");
    std::fs::create_dir_all(&outdir)?;
    index_bam(bam)?;

    let bam_cleaned = remove_hard_clips(bam)?;
    let mut_out = outdir.join(format!("clinvar_{capture}.mut"));
    write_varben_input(&variants, &mut_out)?;
    ensure_bwa_index(&fasta)?;
    let edited_bam = insert_variants(mut_out, bam_cleaned, outdir.clone(), &fasta, mindepth)?;

    write_varben_as_vcf(bam, header, outdir, &fasta)?;

    Ok(edited_bam)
}

/// Use varben (muteditor) to insert a list of variant in a bed files. Require varben, samtools, bwa
/// Require a reference genome and a bwa index
/// Output folder is the varben subfolder of `outdir`
fn insert_variants(
    mut_file: PathBuf,
    bam: PathBuf,
    outdir: PathBuf,
    fasta: &PathBuf,
    mindepth: Option<u32>,
) -> Result<PathBuf, Box<dyn Error>> {
    let outdir_ = outdir.join("varben");
    std::fs::create_dir_all(&outdir_)?;

    let fasta_str = fasta.to_str().ok_or("Invalid fasta path")?;
    let mut_str = mut_file.to_str().ok_or("Invalid mutation file path")?;
    let bam_str = bam.to_str().ok_or("Invalid BAM path")?;
    let outdir_str = outdir_.to_str().ok_or("Invalid output folder")?;

    let output = outdir_.join("edit.sorted.bam");
    // If there's file from older run, override them
    let stale = output.exists()
        && mut_file.metadata().ok().and_then(|m| m.modified().ok())
            > output.metadata().ok().and_then(|m| m.modified().ok());
    if output.exists() && !stale {
        log::debug!("{:?} already exists", output);
    } else {
        let log_path = outdir_.join("varben.log");
        let status = insert_variants_varben(mut_str, bam_str, fasta_str, outdir_str, mindepth, log_path)?;
        if !status.success() {
            return Err(format!("muteditor exited with status {status}").into());
        }
    }
    Ok(output)
}

fn insert_variants_varben(
    mut_str: &str,
    bam_str: &str,
    fasta_str: &str,
    outdir_str: &str,
    mindepth: Option<u32>,
    log_path: PathBuf,
) -> Result<ExitStatus, std::io::Error> {
    let log_file = File::create(&log_path)?;

    let mut cmd = Command::new("muteditor");
    cmd.args([
        "-m", mut_str,
        "-p", "1",
        "-b", bam_str,
        "-r", fasta_str,
        "--aligner", "bwa",
        "--alignerIndex", fasta_str,
        "-o", outdir_str,
    ]);
    if let Some(depth) = mindepth {
        cmd.args(["--mindepth", &depth.to_string()]);
    }
    cmd.stdout(Stdio::from(log_file)).status()
}
/// samtools fastq is the best tool in our testing. It requires read to be sorted by read name beforehand
/// Fastq output will be in the same directory as the BAM and suffixed with _1.fq.gz
pub fn bam_to_fastq(
    bam: PathBuf,
    fq1: PathBuf,
    fq2: PathBuf,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let cmd = format!(
        "samtools sort -n {bam} -@ {threads} | samtools fastq -1 {fq1} -2 {fq2} -0 /dev/null -s /dev/null -n -@ {threads} -",
        bam = bam.display(),
        fq1 = fq1.display(),
        fq2 = fq2.display(),
        threads = nb_threads(),
    );

    let status = Command::new("sh").args(["-c", &cmd]).status()?;

    if !status.success() {
        return Err(format!("bam to fastq conversion exited with {status}").into());
    }
    Ok((fq1, fq2))
}

fn nb_threads() -> usize {
    thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

/// Convert varben tabular data to a simple VCF.
/// Require a fasta as varben only store alternative allel and not refercen
fn write_varben_as_vcf_single(
    mut_file: PathBuf,
    vcf_out: PathBuf,
    mut header: vcf::Header,
    fasta: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    header.sample_names_mut().insert(String::from("TRUTH"));
    header.formats_mut().insert(
        String::from(sample_key::GENOTYPE),
        Map::<Format>::from(sample_key::GENOTYPE),
    );

    let gt_keys: SampleKeys = [String::from(sample_key::GENOTYPE)].into_iter().collect();

    let mut fasta_reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta)?;

    let reader = BufReader::new(File::open(mut_file)?);
    let mut records: Vec<RecordBuf> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 6 {
            continue;
        }
        let chrom = fields[0];
        let pos: usize = fields[1].parse()?;
        let alt = fields[5].trim();

        let region: noodles::core::Region = format!("{}:{}-{}", chrom, pos, pos).parse()?;
        let record = fasta_reader.query(&region)?;
        let seq = record.sequence().as_ref();
        let ref_base = std::str::from_utf8(&seq[..1])?.to_uppercase();

        let samples = SamplesBuf::new(gt_keys.clone(), vec![vec![Some(SampleValue::from("0/1"))]]);
        let buf = RecordBuf::builder()
            .set_reference_sequence_name(chrom)
            .set_variant_start(noodles::core::Position::try_from(pos)?)
            .set_reference_bases(ref_base)
            .set_alternate_bases(AltBasesBuf::from(vec![alt.to_string()]))
            .set_samples(samples)
            .build();

        records.push(buf);
    }

    sort_by_chromosome(&mut records);

    // let header = vcf::Header::default();
    let mut writer = File::create(vcf_out.clone())
        .map(bgzf::io::Writer::new)
        .map(vcf::io::Writer::new)?;
    writer.write_header(&header)?;
    for rec in &records {
        writer.write_variant_record(&header, rec)?;
    }

    log::debug!("Wrote {} variants to {:?}", records.len(), vcf_out);
    Ok(())
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

        write_simuscop_config(
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
        assert!(content.contains("profile = /data/sample.profile"), "missing profile");
        assert!(content.contains("variation = /data/clinvar.simuscop"), "missing variation");
        assert!(content.contains("target = /data/capture.bed"), "missing target");
        assert!(content.contains("name = agilent-col6a1"), "missing name");
        assert!(content.contains("output = /data/simuscop_out"), "missing output");
        assert!(content.contains("layout = PE"), "missing layout");
        assert!(content.contains("coverage = 50"), "missing coverage");
        assert!(content.contains("threads = "), "missing threads");
    }

    #[test]
    fn simuscop_config_coverage_matches_input() {
        let dir = std::env::temp_dir().join("simuscop_coverage_test");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("test.conf");

        write_simuscop_config(
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
