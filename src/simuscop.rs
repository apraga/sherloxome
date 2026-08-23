use crate::resolve_bam;
use crate::run::{run_from_filename, run_to_string};
use crate::silico::{index_vcf, nb_threads, sort_by_chromosome};
use flate2::Compression;
use flate2::write::GzEncoder;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::header::record::value::{Map, map::Format};
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::io::Write as VariantWrite;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::samples::keys::key as sample_key;
use noodles::vcf::variant::record_buf::samples::{
    Keys as SampleKeys, Samples as SamplesBuf, sample::Value as SampleValue,
};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
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

/// Write the variants used to drive simuReads as a bgzipped VCF, giving a ground-truth
/// reference VCF for benchmarking with hap.py (see `benchmark::silico_run_to_happy`).
/// simuReads inserts every requested variant (unlike varben, which can fail to insert some),
/// so the input variant list doubles as the truth set. All variants are heterozygous,
/// matching the "het" zygosity always emitted by `write_variant`.
pub fn write_as_vcf(
    variants: &[RecordBuf],
    header: &vcf::Header,
    vcf_out: &Path,
) -> Result<(), Box<dyn Error>> {
    let mut header = header.clone();
    header.sample_names_mut().insert(String::from("TRUTH"));
    header.formats_mut().insert(
        String::from(sample_key::GENOTYPE),
        Map::<Format>::from(sample_key::GENOTYPE),
    );

    let gt_keys: SampleKeys = [String::from(sample_key::GENOTYPE)].into_iter().collect();

    let mut records: Vec<RecordBuf> = variants
        .iter()
        .map(|v| {
            let mut rec = v.clone();
            let samples =
                SamplesBuf::new(gt_keys.clone(), vec![vec![Some(SampleValue::from("0/1"))]]);
            *rec.samples_mut() = samples;
            rec
        })
        .collect();
    sort_by_chromosome(&mut records);

    let mut writer = File::create(vcf_out)
        .map(bgzf::io::Writer::new)
        .map(vcf::io::Writer::new)?;
    writer.write_header(&header)?;
    for rec in &records {
        writer.write_variant_record(&header, rec)?;
    }
    // Drop (rather than let it fall out of scope at the end of the function) so the
    // BGZF trailer is flushed to disk before tabix reads the file below.
    drop(writer);

    log::debug!("Wrote {} variants to {:?}", records.len(), vcf_out);
    index_vcf(vcf_out)?;
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
    dbsnp_vcf: PathBuf,
) -> Result<(), Box<dyn Error>> {
    let threads = nb_threads();
    let mut w = BufWriter::new(File::create(config_path)?);
    writeln!(w, "ref = {}", fasta.display())?;
    writeln!(w, "profile = {}", profile_dir.display())?;
    writeln!(w, "variation = {}", variation.display())?;
    writeln!(w, "snp = {}", dbsnp_vcf.display())?;
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
    bam: &Option<String>,
    bed: &PathBuf,
    // capture: &str,
    fasta: &PathBuf,
    vcf: Option<&PathBuf>,
    profile: Option<&PathBuf>,
    variants: &Vec<RecordBuf>,
    header: &vcf::Header,
    outdir: &PathBuf,
    coverage: u32,
    snp_file: PathBuf,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let profile_dir = if let Some(p) = profile {
        p.clone()
    } else if let Some(v) = vcf {
        if let Some(bam_file) = bam {
            let bam_path = resolve_bam(bam_file, outdir)?;
            generate_profile(&bam_path, v, fasta, bed, outdir)?
        } else {
            log::error!("[silico.simuscop] requires a BAM file if no profile is set");
            return Err("[silico.simuscop] requires a BAM file if no profile is set".into());
        }
    } else {
        log::error!("[silico.simuscop] requires either `profile` or `vcf`");
        return Err("[silico.simuscop] requires either `profile` or `vcf`".into());
    };

    // Add a false patient id before extracting
    let profile_stem = profile_dir
        .file_stem()
        .ok_or("Invalid profile path")?
        .to_string_lossy();
    let profile2 = PathBuf::from(format!("nopatient_{}", profile_stem));
    let run = run_from_filename(&profile2).ok_or_else(|| {
        format!(
            "profile file {:?} does not follow the SAMPLE_SEQUENCER_CAPTURE_DEPTHx{{_SILICO}} filenaming scheme",
            profile2
        )
    })?;
    let run_str = run_to_string(&run);

    let base = format!("{run_str}_simuscop");

    std::fs::create_dir_all(outdir)?;
    let variation = outdir.join(format!("{run_str}.simuscop"));
    write_input(variants, &variation, &base)?;

    // `benchmark::silico_run_to_happy` always looks up the truth VCF under data/exp_raw,
    // so write it here
    let truth_vcf = PathBuf::from("data/exp_raw").join(format!("{base}.vcf.gz"));
    std::fs::create_dir_all(truth_vcf.parent().unwrap())?;
    write_as_vcf(variants, header, &truth_vcf)?;

    let config_path = outdir.join(format!("{base}.conf"));
    write_config(
        &config_path,
        fasta,
        &profile_dir,
        &variation,
        bed,
        &base,
        outdir,
        coverage,
        snp_file,
    )?;

    run_simu_reads(&config_path, &base, outdir)
}

/// Run simuReads with the given config. Returns (fq1, fq2) output paths.
fn run_simu_reads(
    config_path: &Path,
    name: &str,
    output_dir: &Path,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let raw_fq1 = output_dir.join(format!("{name}_1.fq"));
    let raw_fq2 = output_dir.join(format!("{name}_2.fq"));
    let gz_fq1 = PathBuf::from(format!("{}.gz", raw_fq1.display()));
    let gz_fq2 = PathBuf::from(format!("{}.gz", raw_fq2.display()));

    if gz_fq1.exists() && gz_fq2.exists() {
        log::info!(
            "Skipping BAM editing as output fastq already exists {:?}, {:?}",
            gz_fq1,
            gz_fq2
        );
        return Ok((gz_fq1, gz_fq2));
    }

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

    let fq1 = compress_fastq(&raw_fq1)?;
    let fq2 = compress_fastq(&raw_fq2)?;
    Ok((fq1, fq2))
}

/// Gzip-compress a plain fastq file produced by simuReads, removing the uncompressed
/// original. Returns the path to the resulting `.fq.gz` file.
fn compress_fastq(fq: &Path) -> Result<PathBuf, Box<dyn Error>> {
    let gz_path = PathBuf::from(format!("{}.gz", fq.display()));

    let mut reader = BufReader::new(File::open(fq)?);
    let mut encoder = GzEncoder::new(File::create(&gz_path)?, Compression::default());
    std::io::copy(&mut reader, &mut encoder)?;
    encoder.finish()?;

    std::fs::remove_file(fq)?;
    Ok(gz_path)
}

/// Build a seqToProfile sequencing profile from a normal BAM.
/// The profile directory is derived from the BAM stem and reused on subsequent runs.
fn generate_profile(
    bam_path: &PathBuf,
    vcf: &PathBuf,
    fasta: &PathBuf,
    bed: &PathBuf,
    outdir: &PathBuf,
) -> Result<PathBuf, Box<dyn Error>> {
    let bam_stem = bam_path.file_stem().unwrap().to_string_lossy();
    let profile_file = outdir.join(format!("{bam_stem}.profile"));

    if profile_file.exists() {
        log::debug!("Reusing existing profile: {:?}", profile_file);
        return Ok(profile_file);
    }
    std::fs::create_dir_all(outdir)?;
    log::info!("Generating sequencing profile with seqToProfile");
    let args = [
        "-b",
        bam_path.to_str().ok_or("Invalid BAM path")?,
        "-v",
        vcf.to_str().ok_or("Invalid VCF path")?,
        "-r",
        fasta.to_str().ok_or("Invalid FASTA path")?,
        "-t",
        bed.to_str().ok_or("Invalid BED path")?,
        "-o",
        profile_file.to_str().ok_or("Invalid profile dir path")?,
    ];
    println!("{:?}", args);

    let status = Command::new("seqToProfile").args(args).status()?;

    if !status.success() {
        // Remove the empty directory so the next run retries cleanly
        let _ = std::fs::remove_dir(&profile_file);
        return Err(format!("seqToProfile exited with status {status}").into());
    }

    Ok(profile_file)
}
