use crate::silico::nb_threads;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record::AlternateBases;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Write one SNV in simuscop tab-separated variation format (no header line).
/// Always emits `het` since we sample heterozygous AF from clinvar.
fn write_variant(
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
pub fn write_input(variants: &[RecordBuf], out: &Path, sample: &str) -> Result<(), Box<dyn Error>> {
    let mut w = BufWriter::new(File::create(out)?);
    for v in variants {
        write_variant(&mut w, v, sample)?;
    }
    w.flush()?;
    Ok(())
}

/// Write a simuReads config file
pub fn write_config(
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

/// Generate controls into a in-silico FASTQ.
/// Either `profile` (pre-built directory) or `vcf` (runs seqToProfile) must be provided.
pub fn generate_controls_fastq(
    bam: &PathBuf,
    bed: &PathBuf,
    capture: &str,
    fasta: &PathBuf,
    vcf: Option<&PathBuf>,
    profile: Option<&PathBuf>,
    variants: &Vec<RecordBuf>,
    outdir: &PathBuf,
    coverage: u32,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let profile_dir = if let Some(p) = profile {
        p.clone()
    } else if let Some(v) = vcf {
        ensure_profile(bam, v, fasta, bed, outdir)?
    } else {
        log::error!("[silico.simuscop] requires either `profile` or `vcf`");
        return Err("[silico.simuscop] requires either `profile` or `vcf`".into());
    };

    let variation = outdir.join(format!("clinvar_{capture}.simuscop"));
    write_input(variants, &variation, capture)?;

    let simuscop_outdir = outdir.join(format!("simuscop_{capture}"));
    std::fs::create_dir_all(&simuscop_outdir)?;

    let config_path = outdir.join(format!("simuscop_{capture}.conf"));
    write_config(
        &config_path,
        fasta,
        &profile_dir,
        &variation,
        bed,
        capture,
        &simuscop_outdir,
        coverage,
    )?;

    run_simu_reads(&config_path, capture, &simuscop_outdir)
}

/// Run simuReads with the given config. Returns (fq1, fq2) output paths.
fn run_simu_reads(
    config_path: &Path,
    name: &str,
    output_dir: &Path,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let status = Command::new("simuReads")
        .arg(
            config_path
                .to_str()
                .ok_or("Invalid simuReads config path")?,
        )
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
