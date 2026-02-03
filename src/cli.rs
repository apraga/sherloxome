use crate::fastqbaid2020::{Depth, Kit, Run, Sequencer};
use crate::fastqbaid2020::{available, samplesheet_real};
use crate::giab::Patient;
use clap::{Parser, Subcommand};
use glob::glob;
use itertools::iproduct;
use std::collections::HashSet;
use std::path::PathBuf;
use toml;

use serde::Deserialize;

/// Configuration file for the user definig which GIAB data and which in silico data
#[derive(Deserialize, Debug)]
struct Config {
    real: Option<RealConfig>,
    silico: Option<SilicoConfig>,
}

#[derive(Deserialize, Debug)]
struct SilicoConfig {
    patients: Vec<Patient>,
}

#[derive(Deserialize, Debug)]
struct RealConfig {
    patients: Vec<Patient>,
    sequencers: Vec<Sequencer>,
    kits: Vec<Kit>,
    depths: Vec<Depth>,
}

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate samplesheet for all use cases
    Samplesheet {
        /// Sets a custom config file (TOML)
        #[arg(short, long, value_name = "FILE")]
        config: PathBuf,
    },
    /// Analyse VCF with hap.py
    Analyse {
        #[arg(short, long, value_name = "INPUT_DIR")]
        input: PathBuf,

        #[arg(short, long, value_name = "OUTPUT_DIR")]
        output: PathBuf,
    },
}

/// Do the cartesian product of the runs select by the user
fn candidate_runs(real: &RealConfig) -> HashSet<Run> {
    iproduct!(&real.patients, &real.sequencers, &real.kits, &real.depths)
        .map(|(p, s, k, d)| Run {
            patient: *p,
            sequencer: *s,
            kit: *k,
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
        real.kits.len(),
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

/// Generate a samplesheet for real and in-silico data according to a configuration file
fn generate_samplesheet(config: &Config) {
    if let Some(real) = &config.real {
        let runs = filter_available_runs(&real);
        samplesheet_real(runs);
    } else {
        println!("No real patients...");
    }

    if let Some(silico) = &config.silico {
        println!("Silico patients: {:?}", silico.patients);
    } else {
        println!("No silico patients...");
    }
}

/// Analyse all VCF in a directory.
///
/// Match them to a reference patients, call happy on each VCF and its
/// reference and merge the results.
fn analyse(input_dir: PathBuf, output_dir: PathBuf) {
    // Convert directory to string, not easy
    let pattern = format!("{}/**/*.vcf.gz", input_dir.display().to_string());
    println!("{:?}", pattern);
    for vcf in glob(pattern.as_str()).unwrap() {
        if let Ok(path) = vcf {
            println!("{:?}", path.display())
        }
    }
}

/// Read CLI arguments and call relevant functions
pub fn process_cli() {
    let cli = Cli::parse();
    match &cli.command {
        Some(Commands::Samplesheet { config }) => {
            let content = std::fs::read_to_string(config).unwrap();
            let c: Config = toml::from_str(&content).unwrap();
            generate_samplesheet(&c);
        }
        Some(Commands::Analyse { input, output }) => {
            analyse(input.to_path_buf(), output.to_path_buf());
        }
        None => {}
    }

    // Continued program logic goes here...
}
