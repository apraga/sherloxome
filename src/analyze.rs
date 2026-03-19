use crate::fastqbaid2020::{Run, run_from_filename};
use crate::giab;
use std::collections::HashMap;
// use flate2::read::GzDecoder;
use glob::glob;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::error::Error;
use std::fs::create_dir_all;
use std::path::PathBuf;
use std::process::Command;

///! Analyze all VCF in a directory with happy by comparing to reference VCF

/// From a directory, list all VCFs matching runs information, call happy on each VCF
/// merge the results.
/// If there's an error in a VCF analysis, print a message and skip to the next
pub fn analyze(
    fasta: &PathBuf,
    capture: HashMap<String, String>,
    input_dir: PathBuf,
    output_dir: PathBuf,
) -> Result<(), Box<dyn Error>> {
    let sdf = fasta_to_sdf(fasta);
    create_dir_all(&output_dir)?;
    vcf_to_analyze(input_dir)
        .into_par_iter()
        .filter_map(|(query_vcf, run)| {
            prepare_analysis(&run, &query_vcf, &capture).map(|(query_bed, truth_vcf, truth_bed)| {
                (query_vcf, run, query_bed, truth_vcf, truth_bed)
            })
        })
        .for_each(|(query_vcf, run, query_bed, truth_vcf, truth_bed)| {
            happy_vcfeval(
                &truth_vcf,
                &truth_bed,
                &query_vcf,
                &query_bed,
                &output_dir,
                fasta,
                &sdf,
                run,
            );
        });
    Ok(())
}

/// Prepare files for happy : construct the filepath and check it exists
fn prepare_analysis(
    run: &Run,
    query_vcf: &PathBuf,
    capture: &HashMap<String, String>,
) -> Option<(PathBuf, PathBuf, PathBuf)> {
    let query_bed_ = capture.get(&run.capture.to_string()).or_else(|| {
        eprintln!("No BED file for capture {:?}", run.capture);
        None
    })?;
    let query_bed = PathBuf::from(query_bed_);
    let truth_vcf = PathBuf::from("data/ref").join(giab::vcf_file(&run.patient));
    let truth_bed = PathBuf::from("data/ref").join(giab::bed_file(&run.patient));
    if let Some(missing) = [&query_bed, &query_vcf, &truth_vcf, &truth_bed]
        .iter()
        .find(|p| !p.exists())
    {
        eprintln!("Missing {:?}, skipping analysis", missing);
        return None;
    }
    Some((query_bed, truth_vcf, truth_bed))
}

/// From a directory, list all VCFs where the filename has run information
/// i.e patient, capture, sequencer, depth. The last 2 parameters are not mandatory per se as GIAB reference
/// don't have it but it's useful for statistical analysis later on.
fn vcf_to_analyze(input_dir: PathBuf) -> Vec<(PathBuf, Run)> {
    let pattern = format!("{}/**/*.vcf.gz", input_dir.display());
    match glob(pattern.as_str()) {
        Err(e) => {
            eprintln!("Invalid glob pattern: {}", e);
            vec![]
        }
        Ok(paths) => paths
            .flatten()
            .filter_map(|path| run_from_filename(&path).map(|run| (path, run)))
            .collect(),
    }
}

/// Generate a sdf file from FASTA for happy. SDF file will be stored in the same directory as the FASTA
fn fasta_to_sdf(fasta: &PathBuf) -> PathBuf {
    let sdf = fasta.with_extension("sdf");
    if !sdf.exists() {
        Command::new("rtg")
            .arg("format")
            .arg(fasta)
            .arg("-o")
            .arg(&sdf)
            .status()
            .expect("Failed to execute command");
    } else {
        println!("{} already exists", sdf.display())
    }
    sdf
}

// Run happy with vcfeval (see for example the nix package)
// We parallelize over the number of VCF. rtg is run as single thread for now.
// Warning : --threads must be set, otherwise the default value on some cluster may cause it to fail.
fn happy_vcfeval(
    truth_vcf: &PathBuf,
    truth_bed: &PathBuf,
    query_vcf: &PathBuf,
    query_bed: &PathBuf,
    output_dir: &PathBuf,
    fasta: &PathBuf,
    sdf: &PathBuf,
    run: Run,
) {
    let prefix = format!(
        "{}-{}-{}-{}",
        run.patient, run.capture, run.sequencer, run.depth
    );
    let out = output_dir.join(prefix);
    if out.join(".summary.csv").exists() {
        println!("Summary file already exists (summary)");
        return;
    }

    println!("Processing {:?}", query_vcf);
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
            // "--sample",
            // &run.patient.to_string(),
            "-o",
            out.to_str().unwrap(),
        ])
        .status()
        .expect("Failed to run hap.py");

    if !status.success() {
        eprintln!("hap.py failed for {:?}", query_vcf);
    }
}
