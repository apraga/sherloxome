//! Setup data before the pipeline for either real patients or insilico data:
//! - download all reference FASTQ for real patients and BED in a directory.
//! - generate in silico FASTQ from clinvar data [`crate::silico`] by either
//!   - inserting variants into a real BAM
//!   - generatingfFASTQ directly from a list of variants (simuscop)
//! - generate a samplesheet for sarek.
//!
//! At least one configuration should be set (real patients or insilico data)

use crate::baid2020::{Capture, Sequencer};
use crate::baid2020::{available, real_row};
use crate::download_blocking;
use crate::giab::{Patient, bed_url, tbi_url, vcf_url};
use crate::giab::{bed_filename, tbi_path, vcf_filename};
use crate::resolve_fasta;
use crate::run::Run;
use crate::silico::*;
use itertools::iproduct;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

/// Configuration file for the user definig which GIAB data and which in silico data
#[derive(Deserialize, Debug)]
pub struct Config {
    pub fasta: String,
    real: Option<RealConfig>,
    silico: Option<SilicoConfig>,
    pub capture: HashMap<String, String>,
}

impl Config {
    pub fn validate(&self) -> Result<(), Box<dyn Error>> {
        let mut required: Vec<String> = Vec::new();
        if let Some(real) = &self.real {
            required.extend(real.captures.iter().map(|c| c.to_string()));
        }
        if let Some(silico) = &self.silico {
            required.push(silico.capture.clone());
        }
        let missing: Vec<String> = required
            .into_iter()
            .filter(|name| !self.capture.contains_key(name))
            .collect();
        if !missing.is_empty() {
            return Err(format!(
                "Missing required captures in [capture] section: {}",
                missing.join(", ")
            )
            .into());
        }
        Ok(())
    }
}

#[derive(Deserialize, Debug)]
pub struct RealConfig {
    pub patients: Vec<Patient>,
    pub sequencers: Vec<Sequencer>,
    pub captures: Vec<Capture>,
    pub depths: Vec<u32>,
}

pub struct SamplesheetRow {
    pub patient: String,
    pub sample: String,
    pub lane: u8,
    pub fastq_1: String,
    pub fastq_2: String,
}

/// Write all rows to `samplesheet.csv` in the current directory.
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
        sample: format!("{}", *p),
        sequencer: (*s).to_string(),
        capture: (*c).to_string(),
        depth: *d,
        silico: None,
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
pub fn download_giab_run(patient: &Patient) {
    let out_dir = PathBuf::from("data/ref");
    std::fs::create_dir_all(&out_dir).expect("Could not create directory");
    let vcf = out_dir.join(vcf_filename(patient));
    let bed = out_dir.join(bed_filename(patient));
    download_blocking(&vcf_url(patient), &vcf);
    download_blocking(&tbi_url(patient), &tbi_path(patient));
    download_blocking(&bed_url(patient), &bed);
}

fn download_giab_runs(runs: &HashSet<Run>) {
    let patients: HashSet<Patient> = runs
        .iter()
        .map(|r| {
            r.sample
                .parse::<Patient>()
                .expect("Failed to convert sample to patient")
        })
        .collect();
    for p in patients {
        download_giab_run(&p);
    }
}

/// Main entry point for the `setup` subcommand.
///
/// Depending on what is present in `conf`, downloads GIAB data and/or generates
/// in silico controls, then writes `samplesheet.csv`.
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
        let capture_str = &silico.capture;
        let bed = conf
            .capture
            .get(capture_str)
            .ok_or("no BED for silico capture")?;
        let fasta = resolve_fasta(&conf.fasta)?;
        let rows_silico = generate_controls(silico, PathBuf::from(bed), capture_str, fasta)?;
        rows.extend(rows_silico);
    }
    write_samplesheet(rows)?;
    Ok(())
}

/// Build a samplesheet row for an in silico sample.
pub fn silico_row(silico_type: &str, fq1: PathBuf, fq2: PathBuf) -> SamplesheetRow {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::baid2020::{Capture, Sequencer};
    use crate::giab::Patient;

    fn make_real(captures: Vec<Capture>) -> RealConfig {
        RealConfig {
            patients: vec![Patient::HG002],
            sequencers: vec![Sequencer::Novaseq],
            captures,
            depths: vec![50],
        }
    }

    fn make_silico(capture: &str) -> SilicoConfig {
        SilicoConfig {
            capture: capture.to_string(),
            bam_file: String::new(),
            clinvar: None,
            nb_variants: None,
            outdir: None,
            simuscop: None,
            varben: None,
        }
    }

    fn config(
        real: Option<RealConfig>,
        silico: Option<SilicoConfig>,
        captures: &[(&str, &str)],
    ) -> Config {
        Config {
            fasta: String::new(),
            real,
            silico,
            capture: captures
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect(),
        }
    }

    #[test]
    fn accepts_empty_config() {
        assert!(config(None, None, &[]).validate().is_ok());
    }

    #[test]
    fn rejects_real_with_missing_capture() {
        let conf = config(
            Some(make_real(vec![Capture::Agilent, Capture::Idt])),
            None,
            &[("agilent", "agilent.bed")],
        );
        let msg = conf.validate().unwrap_err().to_string();
        assert!(msg.contains("idt"), "expected idt in error: {msg}");
    }

    #[test]
    fn accepts_real_with_all_captures_present() {
        let conf = config(
            Some(make_real(vec![Capture::Agilent, Capture::Idt])),
            None,
            &[("agilent", "agilent.bed"), ("idt", "idt.bed")],
        );
        assert!(conf.validate().is_ok());
    }

    #[test]
    fn rejects_silico_with_missing_capture() {
        let conf = config(None, Some(make_silico("agilent-col6a1")), &[]);
        let msg = conf.validate().unwrap_err().to_string();
        assert!(
            msg.contains("agilent-col6a1"),
            "expected agilent-col6a1 in error: {msg}"
        );
    }

    #[test]
    fn accepts_silico_with_custom_capture() {
        let conf = config(
            None,
            Some(make_silico("agilent-col6a1")),
            &[("agilent-col6a1", "col6a1.bed")],
        );
        assert!(conf.validate().is_ok());
    }
}
