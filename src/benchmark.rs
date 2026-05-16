//! Benchmark variant-called VCFs against GIAB truth sets using hap.py.
//!
//! Run metadata (patient, sequencer, capture kit, depth) is extracted from VCF filenames,
//! so files must be named using the conventions established in [`crate::fastqbaid2020`].
use crate::baid2020;
use crate::check_deps;
use crate::giab;
use crate::ref_dir;
use crate::resolve_fasta;
use crate::run::{Run, run_from_filename, run_to_string};
use crate::setup::Config; //, RealConfig, SilicoConfig, filter_available_runs};
use glob::glob;
use polars::lazy::prelude::*;
use polars::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
// use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::create_dir_all;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Files needed for happy
struct HappyRunSetup {
    prefix: String,
    query_vcf: PathBuf,
    query_bed: PathBuf,
    truth_vcf: PathBuf,
    truth_bed: PathBuf,
}

/// Analyze VCFs derived from config against their truth VCF with hap.py, then merge results.
/// Real patient VCFs use GIAB as truth; the silico VCF uses the `_success.vcf.gz` from setup.
pub fn analyze(
    conf: &Config,
    // fasta: &str,
    // capture: HashMap<String, String>,
    input_dir: PathBuf,
    output_dir: PathBuf,
) -> Result<(), Box<dyn Error>> {
    check_deps(&["rtg"]);
    create_dir_all(&output_dir)?;
    let fasta = resolve_fasta(&conf.fasta)?;
    let sdf = fasta_to_sdf(&fasta)?;
    let runs = runs_to_analyze(&input_dir)?;
    let queue: Vec<HappyRunSetup> = runs
        .into_iter()
        .filter_map(|(run, vcf)| run_to_happy(run, vcf))
        .collect();

    download_missing(&queue);

    queue
        .into_par_iter()
        .filter_map(|item| {
            validate_files(&[
                &item.query_vcf,
                &item.query_bed,
                &item.truth_vcf,
                &item.truth_bed,
            ])
            .map(|_| item)
        })
        .for_each(|item| {
            happy_vcfeval(
                &item.truth_vcf,
                &item.truth_bed,
                &item.query_vcf,
                &item.query_bed,
                &output_dir,
                &fasta,
                &sdf,
                &item.prefix,
            );
        });

    merge_summaries(output_dir)?;
    Ok(())
}

/// For all VCFs in a directory matching the filenaming scheme in `crate::run`, find a corresponding `Run`
fn runs_to_analyze(input_dir: &Path) -> Result<Vec<(Run, PathBuf)>, Box<dyn Error>> {
    let pattern = format!("{}/**/*.vcf.gz", input_dir.display());
    let candidates: Vec<PathBuf> = glob(&pattern)?.into_iter().filter_map(|p| p.ok()).collect();
    log::debug!("Found {} compressed vcf", candidates.len());
    if candidates.len() == 0 {
        return Err(format!("No .vcf.gz found in {} !", input_dir.display()).into());
    }
    let runs: Vec<(Run, PathBuf)> = candidates
        .into_iter()
        .filter_map(|p| {
            let filename = p.file_name().map(PathBuf::from)?;
            let run = run_from_filename(&filename)?;
            Some((run, p))
        })
        .collect();
    if runs.len() == 0 {
        return Err(format!("No VCF matches the filenaming scheme !").into());
    }
    log::info!("{} runs to analyze", runs.len());
    Ok(runs)
}

/// For a run, find reference VCF and capture for happy
fn run_to_happy(run: Run, query_vcf: PathBuf) -> Option<HappyRunSetup> {
    match run.silico {
        Some(_) => silico_run_to_happy(run, query_vcf).ok(),
        None => real_run_to_happy(run, query_vcf)
            .map_err(|e| log::warn!("{e}"))
            .ok(),
    }
}

/// For a GIB run, find reference data
fn real_run_to_happy(run: Run, query_vcf: PathBuf) -> Result<HappyRunSetup, Box<dyn Error>> {
    let patient = run
        .sample
        .parse::<giab::Patient>()
        .map_err(|_| format!("unknown GIAB patient: {}", run.sample))?;
    let capture = run
        .capture
        .parse::<baid2020::Capture>()
        .map_err(|_| format!("unknown capture: {}", run.capture))?;
    let truth_vcf = giab::vcf_path(&patient);
    let truth_bed = giab::bed_path(&patient);
    let query_bed = baid2020::capture_path(capture);
    Ok(HappyRunSetup {
        prefix: run_to_string(&run),
        query_vcf,
        query_bed,
        truth_vcf,
        truth_bed,
    })
}

/// TODO support fastq generation with simuscop
/// For silico run, the reference vcf is simply SAMPLE_success.vcf.gz
/// There is no truth_bed, we simple use the capture
fn silico_run_to_happy(run: Run, query_vcf: PathBuf) -> Result<HappyRunSetup, Box<dyn Error>> {
    let prefix = run_to_string(&run);
    let truth_vcf = ref_dir().join(run.sample).join("_success.vcf.gz");
    let capture = run
        .capture
        .parse::<baid2020::Capture>()
        .map_err(|_| format!("unknown capture: {}", run.capture))?;
    let query_bed = baid2020::capture_path(capture);
    let truth_bed = query_bed.clone();
    Ok(HappyRunSetup {
        prefix,
        query_vcf,
        query_bed,
        truth_vcf,
        truth_bed,
    })
}

/// Extract the filename stem from a BAM URL or local path (strips directory and `.bam`).
fn bam_stem(bam: &str) -> String {
    let filename = if bam.starts_with("http://") || bam.starts_with("https://") {
        bam.split('/').last().unwrap_or(bam).to_string()
    } else {
        let path = PathBuf::from(bam);
        path.file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(bam)
            .to_string()
    };
    filename.trim_end_matches(".bam").to_string()
}

/// Merge all hap.py summary CSVs in output_dir into a single merged.csv with run metadata columns.
fn merge_summaries(output_dir: PathBuf) -> Result<(), Box<dyn Error>> {
    let pattern = format!("{}/*.summary.csv", output_dir.display());

    let dfs: Vec<LazyFrame> = glob(&pattern)?
        .flatten()
        .map(|path| -> Result<LazyFrame, Box<dyn Error>> {
            let run =
                run_from_filename(&path).ok_or("Failed to extract run from summary filename")?;
            Ok(CsvReadOptions::default()
                .with_has_header(true)
                .try_into_reader_with_file_path(Some(path))?
                .finish()?
                .lazy()
                .with_column(lit(run.sample).alias("patient"))
                .with_column(lit(run.capture.to_string()).alias("capture"))
                .with_column(lit(run.sequencer.to_string()).alias("sequencer"))
                .with_column(lit(run.depth.to_string()).alias("depth")))
        })
        .collect::<Result<_, _>>()?;

    let mut merged = concat(dfs, UnionArgs::default())?.collect()?;
    let out_file = output_dir.join("merged.csv");
    let mut file = std::fs::File::create(out_file.clone())?;
    CsvWriter::new(&mut file).finish(&mut merged)?;
    println!("Merged results in {:?}", out_file);
    Ok(())
}

/// Generate a SDF index from FASTA for hap.py (stored alongside the FASTA file).
fn fasta_to_sdf(fasta: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    let sdf = fasta.with_extension("sdf");
    if !sdf.join("mainIndex").exists() {
        Command::new("rtg")
            .arg("format")
            .arg(fasta)
            .arg("-o")
            .arg(&sdf)
            .status()
            .map_err(|e| format!("Failed to execute rtg format: {e}"))?;
    } else {
        println!("{} already exists", sdf.display())
    }
    Ok(sdf)
}

fn happy_vcfeval(
    truth_vcf: &PathBuf,
    truth_bed: &PathBuf,
    query_vcf: &PathBuf,
    query_bed: &PathBuf,
    output_dir: &PathBuf,
    fasta: &PathBuf,
    sdf: &PathBuf,
    prefix: &str,
) {
    let out = output_dir.join(prefix);
    if out.with_extension("summary.csv").exists() {
        println!("Summary file already exists, skipping");
        return;
    }

    // println!("Processing {:?}", query_vcf);
    let status = Command::new("hap.py")
        .args([
            truth_vcf.to_str().unwrap(),
            query_vcf.to_str().unwrap(),
            "--false-positives",
            truth_bed.to_str().unwrap(),
            "--target-regions",
            query_bed.to_str().unwrap(),
            "--reference",
            fasta.to_str().unwrap(),
            "--engine=vcfeval",
            "--engine-vcfeval-template",
            sdf.to_str().unwrap(),
            "--threads",
            "1",
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run hap.py");

    if !status.success() {
        eprintln!("hap.py failed for {:?}", query_vcf);
    }
}

/// Download missing GIAB truth VCF (.vcf.gz + .tbi) and BED for every real run in the queue.
/// Silico runs are skipped — their truth VCF is produced by the setup step.
fn download_missing(items: &[HappyRunSetup]) {
    log::debug!("Downloading GIAB missing VCF/BEDs if needed");
    for item in items {
        let Some(patient) = giab::patient_from_filename(&item.truth_vcf) else {
            continue;
        };
        crate::setup::download_giab_run(&patient);
    }
}

fn validate_files(paths: &[&PathBuf]) -> Option<()> {
    if let Some(missing) = paths.iter().find(|p| !p.exists()) {
        log::warn!("Missing {}, skipping analysis", missing.display());
        return None;
    }
    Some(())
}
