//! ## Varben insilico controls
//!
//! Generate control variants by sampling clinvar and insert them into a real BAM

use crate::silico::sort_by_chromosome;
use noodles::bgzf;
use noodles::fasta;
use noodles::vcf;
use noodles::vcf::header::record::value::{Map, map::Format};
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::io::Write as VariantWrite;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::samples::keys::key as sample_key;
use noodles::vcf::variant::record_buf::AlternateBases as AltBasesBuf;
use noodles::vcf::variant::record_buf::samples::{
    Keys as SampleKeys, Samples as SamplesBuf, sample::Value as SampleValue,
};
use rand::RngExt;
use std::error::Error;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;
use std::process::ExitStatus;
use std::process::{Command, Stdio};

/// Write a mutation file for varben using sampled clinvar data
pub fn write_input(variants: &Vec<RecordBuf>, mut_out: &PathBuf) -> Result<(), Box<dyn Error>> {
    // Prepare mutation file
    let mut varben_w = BufWriter::new(File::create(&mut_out)?);
    writeln!(varben_w, "#chrom\tstart\tend\tAF\ttype\talt")?;

    for v in variants {
        write_output(&mut varben_w, v)?;
    }
    varben_w.flush()?;
    Ok(())
}

/// Write varben output as VCF successful insertions only, failures are not properly formatted,
pub fn write_as_vcf(
    bam: &PathBuf,
    header: vcf::Header,
    outdir: PathBuf,
    fasta: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let varben_dir = outdir.join("varben");
    let bam_stem = bam.file_stem().unwrap().to_string_lossy();

    let success_list = varben_dir.join("success_list.txt");
    let success_vcf = outdir.join(format!("{bam_stem}_varben.vcf.gz"));
    write_as_vcf_single(success_list, success_vcf, header, fasta)
}

/// Convert varben tabular data to a simple VCF.
/// Require a fasta as varben only store alternative allel and not refercen
fn write_as_vcf_single(
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

/// Write VCF record to varben format. Assume chromosoes do NOT have a chr prefix
/// Output generate heterozygous variant (AF betwen 0.4 and 0.6)
fn write_output(w: &mut impl Write, record: &RecordBuf) -> Result<(), Box<dyn Error>> {
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

/// Use varben (muteditor) to insert a list of variant in a bed files. Require varben, samtools, bwa
/// Require a reference genome and a bwa index
/// Output folder is the varben subfolder of `outdir`
pub fn insert_variants(
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
        let status =
            run_insert_variants(mut_str, bam_str, fasta_str, outdir_str, mindepth, log_path)?;
        if !status.success() {
            return Err(format!("muteditor exited with status {status}").into());
        }
    }
    Ok(output)
}

fn run_insert_variants(
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
        "-m",
        mut_str,
        "-p",
        "1",
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
    ]);
    if let Some(depth) = mindepth {
        cmd.args(["--mindepth", &depth.to_string()]);
    }
    cmd.stdout(Stdio::from(log_file)).status()
}
