//! Setup data before the pipeline.
//!  Each step is optional:
//! - download all reference FASTQ for real patients and BED in a directory.
//! - generate in silico FASTQ from clinvar data
//!   - insert variants into a real BAM
//!   - generate FASTQ directly from a list of variants (simuscop)
//! And generate a samplesheet for sarek

use crate::download_blocking;
use crate::fastqbaid2020::{Capture, Depth, Run, Sequencer};
use crate::fastqbaid2020::{available, real_row};
use crate::giab::{Patient, bed_url, vcf_url};
use crate::giab::{bed_file, vcf_file};
use crate::silico::*;
use itertools::iproduct;
// use polars::io::resolve_homedir;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Configuration file for the user definig which GIAB data and which in silico data
#[derive(Deserialize, Debug)]
pub struct Config {
    pub fasta: PathBuf,
    real: Option<RealConfig>,
    silico: Option<SilicoConfig>,
    pub capture: HashMap<String, String>,
}

#[derive(Deserialize, Debug)]
struct SilicoConfig {
    capture: Capture,
    bam: Option<String>,
    /// Local ClinVar VCF path; if absent the file is downloaded from NCBI.
    clinvar: Option<PathBuf>,
    /// Number of clinvar variants to sample to insert in the BAM
    nb_variants: Option<u32>,
    /// Output directory for intermediate files; defaults to "data/exp_raw".
    outdir: Option<PathBuf>,
}

#[derive(Deserialize, Debug)]
struct RealConfig {
    patients: Vec<Patient>,
    sequencers: Vec<Sequencer>,
    captures: Vec<Capture>,
    depths: Vec<Depth>,
}

pub struct SamplesheetRow {
    pub patient: String,
    pub sample: String,
    pub lane: u8,
    pub fastq_1: String,
    pub fastq_2: String,
}

fn write_samplesheet(rows: Vec<SamplesheetRow>) -> Result<(), Box<dyn Error>> {
    let mut file = File::create("samplesheet.csv")?;
    writeln!(file, "patient,sample,lane,fastq_1,fastq_2")?;
    for r in rows {
        writeln!(
            file,
            "{},{},{},{},{}",
            r.patient, r.sample, r.lane, r.fastq_1, r.fastq_2
        )?;
    }
    Ok(())
}
/// Do the cartesian product of the runs select by the user
fn candidate_runs(real: &RealConfig) -> HashSet<Run> {
    iproduct!(
        &real.patients,
        &real.sequencers,
        &real.captures,
        &real.depths
    )
    .map(|(p, s, c, d)| Run {
        patient: *p,
        sequencer: *s,
        capture: *c,
        depth: *d,
    })
    .collect::<HashSet<Run>>()
}

/// For all fastq corresponding to real patients, only keep the one where we actually have data
fn filter_available_runs(real: &RealConfig) -> HashSet<Run> {
    let runs_candidates = candidate_runs(&real);

    println!(
        "You asked for {} runs ({:?} patients x {} kits x {} depths x {} sequencers)",
        runs_candidates.len(),
        real.patients.len(),
        real.captures.len(),
        real.depths.len(),
        real.sequencers.len(),
    );
    let runs_available = available();
    // Filter only available runs
    let runs: HashSet<Run> = runs_candidates
        .intersection(&runs_available)
        .cloned()
        .collect();
    println!("Only {:} are available", runs.len());
    runs
}

/// Download reference VCF and BED for the patients present in the given runs.
fn download_giab_runs(runs: &HashSet<Run>) {
    let out_dir = PathBuf::from("data/ref");
    std::fs::create_dir_all(&out_dir).expect("Could not create directory");
    let patients: HashSet<Patient> = runs.iter().map(|r| r.patient).collect();
    for p in patients {
        println!("{:?}", p);
        let mut vcf = out_dir.clone().join(vcf_file(&p));
        let bed = out_dir.clone().join(bed_file(&p));
        download_blocking(&vcf_url(&p), &vcf);
        vcf.set_extension("gz.tbi");
        download_blocking(&vcf_url(&p), &vcf);
        download_blocking(&bed_url(&p), &bed);
    }
}

/// Generate controls from clinvar data and either a BAM file (real patient) or 100% in silico
/// Returns a list of samplesheet rows for writing
fn generate_controls(
    silico: &SilicoConfig,
    bed: PathBuf,
    capture: &str,
    fasta: PathBuf,
) -> Result<Vec<SamplesheetRow>, Box<dyn Error>> {
    check_controls_deps();

    let outdir = silico
        .outdir
        .clone()
        .unwrap_or_else(|| PathBuf::from("data/exp_raw"));

    let mut rows: Vec<SamplesheetRow> = Vec::new();
    if let Some(bam) = &silico.bam {
        let (fq1, fq2) = generate_controls_bam(
            bam,
            bed,
            capture,
            fasta,
            silico.clinvar.clone(),
            silico.nb_variants,
            outdir,
        )?;
        rows.push(silico_row("varben", fq1, fq2));
    } else {
        return Err("FASTQ generation not implemented".into());
    }
    Ok(rows)
}

/// Generate controls from a BAM file and returns 2 fastq
fn generate_controls_bam(
    bam: &str,
    bed: PathBuf,
    capture: &str,
    fasta: PathBuf,
    clinvar: Option<PathBuf>,
    nb_variants: Option<u32>,
    outdir: PathBuf,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    std::fs::create_dir_all(&outdir)?;

    let bam_path = resolve_bam(bam, &outdir)?;
    let bam_stem = bam_path.file_stem().unwrap().to_string_lossy();
    let bam_parent = bam_path.parent().unwrap_or(Path::new("."));
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
        let new_bam = edit_bam(
            &bam_path,
            &bed,
            capture,
            fasta,
            clinvar,
            nb_variants,
            outdir,
        )?;
        bam_to_fastq(new_bam, fq1, fq2)
    }
}

pub fn setup(conf: Config) -> Result<(), Box<dyn Error>> {
    let mut rows: Vec<SamplesheetRow> = Vec::new();
    if conf.real.is_none() && conf.silico.is_none() {
        return Err("Nothing to do: config must have at least one of [real] or [silico]".into());
    }

    if let Some(real) = &conf.real {
        log::debug!("Real patients selected");
        let runs = filter_available_runs(real);
        download_giab_runs(&runs);
        rows.extend(runs.iter().map(real_row));
    }

    if let Some(silico) = &conf.silico {
        log::debug!("Silico patients selected");
        let capture_str = &silico.capture.to_string();
        let bed = conf
            .capture
            .get(capture_str)
            .ok_or("no BED for silico capture")?;
        let fasta = resolve_fasta(&conf.fasta);
        let rows_silico = generate_controls(silico, PathBuf::from(bed), capture_str, fasta)?;
        rows.extend(rows_silico);
    }
    write_samplesheet(rows)?;
    Ok(())
}

fn silico_row(silico_type: &str, fq1: PathBuf, fq2: PathBuf) -> SamplesheetRow {
    let fq1_str = fq1
        .clone()
        .into_os_string()
        .into_string()
        .expect("Unable to confvert fastq 1 to string");
    let fq2_str = fq2
        .into_os_string()
        .into_string()
        .expect("Unable to confvert fastq 2 to string");
    let sample = fq1
        .file_stem()
        .unwrap()
        .to_os_string()
        .into_string()
        .unwrap()
        .replace("_1.fq", "");

    return SamplesheetRow {
        patient: format!("silico-{}", silico_type),
        lane: 1,
        sample: sample,
        fastq_1: fq1_str,
        fastq_2: fq2_str,
    };
}
