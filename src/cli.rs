use crate::fastqbaid2020::{Depth, Kit, Run, Sequencer};
use crate::fastqbaid2020::{available, samplesheet_real};
use crate::giab::Patient;
use clap::{Parser, Subcommand};
use itertools::iproduct;
use std::collections::HashSet;
use std::path::PathBuf;
use toml;

use serde::Deserialize;

/// Configuration file for the user definig which GIAB data and which in silico data
#[derive(Deserialize, Debug)]
pub struct Config {
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
    /// Sets a custom config file (TOML)
    #[arg(short, long, value_name = "FILE")]
    config: PathBuf,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate samplesheet for all use cases
    Samplesheet {},
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

/// Read CLI arguments and call relevant functions
pub fn process_cli() {
    let cli = Cli::parse();
    let content = std::fs::read_to_string(cli.config).unwrap();

    let config: Config = toml::from_str(&content).unwrap();
    match &cli.command {
        Some(Commands::Samplesheet {}) => {
            generate_samplesheet(&config);
        }
        None => {}
    }

    // Continued program logic goes here...
}
