//! # Insilico Controls
//! Generate control variants by sampling clinvar data

use crate::download_blocking;
use log;
use noodles::bed;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::info::field::Value;
use rand::RngExt;
use rand::prelude::IteratorRandom;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitStatus;
use std::process::{Command, Stdio};
use std::thread;
use url::Url;
use which::which;

struct Snv {
    chrom: String, // chr-prefixed
    pos: u64,      // 1-based
    alt: String,
}

pub fn check_controls_deps() {
    let tools = ["bwa", "samtools", "muteditor"];

    log::debug!("Checking dependencies");
    for tool in tools {
        match which(tool) {
            Ok(path) => log::debug!("{tool} ✓  ({})", path.display()),
            Err(_) => panic!("{tool} ✗  not found in PATH"),
        }
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
/// Assubme bed is stored
fn load_bed(bed: &PathBuf) -> Result<BedIndex, Box<dyn Error>> {
    let mut reader = File::open(bed)
        .map(BufReader::new)
        .map(bed::io::Reader::<3, _>::new)?;

    let mut index = BedIndex::new();
    let mut record = bed::Record::default();

    while reader.read_record(&mut record)? != 0 {
        let chrom = std::str::from_utf8(record.reference_sequence_name())?
            .trim_start_matches("chr")
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

// Generate AF betwen 0.4 and 0.6
fn write_varben(
    w: &mut impl Write,
    snv: &Snv,
    rng: &mut impl RngExt,
) -> Result<(), Box<dyn Error>> {
    let af: f64 = rng.random_range(0.4..0.6);
    writeln!(
        w,
        "{}\t{}\t{}\t{}\tsnv\t{}",
        snv.chrom, snv.pos, snv.pos, af, snv.alt
    )?;
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

fn sample_clinvar_variants(
    reader: &mut vcf::io::Reader<bgzf::io::Reader<File>>,
    header: &vcf::Header,
    capture: &BedIndex,
    spacing: u32,
    n: usize,
    rng: &mut impl RngExt,
) -> Vec<(vcf::Record, Snv)> {
    let mut last_pos: HashMap<String, u64> = HashMap::new();
    let mut selected: Vec<(vcf::Record, Snv)> = reader
        .records()
        .filter_map(|r| r.ok())
        .filter_map(|record| {
            keep_variant(&record, &header, &capture, &last_pos, spacing).map(|snv| {
                last_pos.insert(snv.chrom.clone(), snv.pos);
                (record, snv)
            })
        })
        .sample(rng, n);
    selected.sort_by(|(_, a), (_, b)| {
        chrom_order(&a.chrom)
            .cmp(&chrom_order(&b.chrom))
            .then(a.pos.cmp(&b.pos))
    });

    log::debug!("Selected {} variants", selected.len());
    selected
}

/// Select n clinvar pathogenic SNV inside the capture kit 50bp apart.
/// Output is both a mutation file for varben and a VCF (hardocoded into data/exp_raw/clinvar_$CAPTURE)
/// The most efficient way is to parse clinvar once to get variant in the capture kit and 50pb apart.
/// In a second pass, sample randomly n of theme
pub fn sample_clinvar(
    bed: PathBuf,
    spacing: u32,
    n: usize,
    vcf_out: PathBuf,
    mut_out: PathBuf,
) -> Result<(), Box<dyn Error>> {
    let vcf = download_clinvar();

    // Prepare VCF file
    // Load BED into interval tree per chrom
    let capture = load_bed(&bed)?;
    let mut reader = File::open(vcf)
        .map(bgzf::io::Reader::new)
        .map(vcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut rng = rand::rng();
    let variants = sample_clinvar_variants(&mut reader, &header, &capture, spacing, n, &mut rng);

    write_sampled_clinvar_vcf(&variants, &header, &vcf_out)?;
    write_sampled_clinvar_mut(&variants, &mut_out, &mut rng)?;
    log::debug!("Wrote {:?} and {:?}", vcf_out, mut_out);
    Ok(())
}

fn write_sampled_clinvar_vcf(
    variants: &Vec<(vcf::Record, Snv)>,
    header: &vcf::Header,
    vcf_out: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let mut writer = File::create(&vcf_out)
        .map(bgzf::io::Writer::new)
        .map(vcf::io::Writer::new)?;
    writer.write_header(&header)?;

    for (record, _) in variants {
        writer.write_record(&header, record)?;
    }
    Ok(())
}

fn write_sampled_clinvar_mut(
    variants: &Vec<(vcf::Record, Snv)>,
    mut_out: &PathBuf,
    rng: &mut impl RngExt,
) -> Result<(), Box<dyn Error>> {
    // Prepare mutation file
    let mut varben_w = BufWriter::new(File::create(&mut_out)?);
    writeln!(varben_w, "#chrom\tstart\tend\tAF\ttype\talt")?;

    for (_, snv) in variants {
        write_varben(&mut varben_w, snv, rng)?;
    }
    varben_w.flush()?;
    Ok(())
}
/// Do it for a single record
/// Input VCF may use chr prefix, output will have `chr` prefix.
fn keep_variant(
    record: &vcf::Record,
    header: &vcf::Header,
    capture: &BedIndex,
    last_pos: &HashMap<String, u64>,
    spacing: u32,
) -> Option<Snv> {
    let chrom_raw = record.reference_sequence_name().to_string();
    // Remove chr prefix if needed
    let chrom_nb = chrom_raw.strip_prefix("chr").unwrap_or(&chrom_raw);
    let chrom_chr = format!("chr{}", chrom_nb);

    let pos = usize::from(record.variant_start()?.ok()?) as u64;
    let last = last_pos.get(chrom_nb).copied().unwrap_or(0);
    let alt = record
        .alternate_bases()
        .iter()
        .filter_map(|a| a.ok())
        .next()
        .unwrap_or(".")
        .to_string();
    (in_capture(capture, &chrom_nb, pos)
        && is_not_benign(&record.info(), header)
        && is_snv(record)
        && pos != last
        && pos.abs_diff(last) >= spacing.into())
    .then_some(Snv {
        chrom: chrom_chr,
        pos,
        alt,
    })
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

/// Download a bam file if it's an url, otherwise check the file exist
pub fn resolve_bam(bam: &str, outdir: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    if bam.starts_with("http://") || bam.starts_with("https://") {
        download_bam(bam, outdir)
    } else {
        let path = PathBuf::from(bam);
        if !path.exists() {
            return Err(format!("BAM file not found: {}", path.display()).into());
        }
        Ok(path)
    }
}

fn download_bam(url: &str, outdir: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    let parsed = Url::parse(&url)?;
    log::info!("Downloading bam file {}", url);
    let filename = parsed
        .path_segments()
        .and_then(|s| s.last())
        .filter(|s| !s.is_empty())
        .ok_or("could not extract filename from URL")?;

    let output = outdir.join(filename);
    download_blocking(&url, &output);
    Ok(output)
}
/// Generate controls from clinvar variant into a BAM file
/// 1. Sample clinvar randomly
/// 2. Generate mutation files from those clinvar variant for varben
/// 2. Apply varben on a single bam file
pub fn edit_bam(
    bam: &PathBuf,
    bed: &PathBuf,
    capture: &str,
    fasta: PathBuf,
) -> Result<PathBuf, Box<dyn Error>> {
    let outdir = PathBuf::from("data/exp_raw");
    std::fs::create_dir_all(&outdir)?;

    index_bam(&bam)?;
    let bam_cleaned = remove_hard_clips(&bam)?;

    // sample 1000 clinvar variant 50bp apart
    let vcf_out = outdir.join(format!("clinvar_{capture}.vcf.gz"));
    let mut_out = outdir.join(format!("clinvar_{capture}.mut"));
    sample_clinvar(PathBuf::from(bed), 50, 1000, vcf_out, mut_out.clone())?;
    insert_variants(mut_out, bam_cleaned, fasta, outdir)
}

/// Use varben (muteditor) to insert a list of variant in a bed files. Require varben, samtools, bwa
/// Require a reference genome and a bwa index
/// Output folder is the varben subfolder of `outdir`
fn insert_variants(
    mut_file: PathBuf,
    bam: PathBuf,
    fasta: PathBuf,
    outdir: PathBuf,
) -> Result<PathBuf, Box<dyn Error>> {
    let outdir_ = outdir.join("varben");
    std::fs::create_dir_all(&outdir_)?;

    let fasta_str = fasta.to_str().ok_or("Invalid fasta path")?;
    let mut_str = mut_file.to_str().ok_or("Invalid mutation file path")?;
    let bam_str = bam.to_str().ok_or("Invalid BAM path")?;
    let outdir_str = outdir_.to_str().ok_or("Invalid output folder")?;

    let output = outdir_.join("edit.sorted.bam");
    if output.exists() {
        log::debug!("{:?} already exists", output);
    } else {
        let log_path = outdir_.join("varben.log");
        let status = insert_variants_varben(mut_str, bam_str, fasta_str, outdir_str, log_path)?;
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
    log_path: PathBuf,
) -> Result<ExitStatus, std::io::Error> {
    let log_file = File::create(&log_path)?;

    Command::new("muteditor")
        .args([
            "-m",
            mut_str,
            "-p",
            &nb_threads().to_string(),
            "-b",
            bam_str,
            "-r",
            fasta_str,
            "--aligner",
            "bwa",
            "--alignerIndex",
            fasta_str,
            "-o",
            outdir_str,
        ])
        .stdout(Stdio::from(log_file))
        .status()
}
/// samtools fastq is the best tool in our testing. It requires read to be sorted by read name beforehand
/// Fastq output will be in the same directory as the BAM and suffixed with _1.fq.gz
pub fn bam_to_fastq(bam: PathBuf, old_bam: PathBuf) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let stem = old_bam.file_stem().unwrap().to_string_lossy();
    let parent = old_bam.parent().unwrap_or(Path::new("."));
    let fq1 = parent.join(format!("{stem}_1.fq.gz"));
    let fq2 = parent.join(format!("{stem}_2.fq.gz"));
    if fq1.exists() && fq2.exists() {
        log::debug!("Output fastq already exists");
        Ok((fq1, fq2))
    } else {
        run_bam_to_fastq(bam, fq1, fq2)
    }
}
pub fn run_bam_to_fastq(
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
