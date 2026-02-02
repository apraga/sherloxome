use clap::{Parser, Subcommand};
use itertools::iproduct;
use sherloxome::fastqbaid2020::{Depth, Kit, Run, Sequencer, available, write_samplesheet};
use sherloxome::giab::Patient;
use std::collections::HashSet;
use std::path::PathBuf;
use toml;

use serde::Deserialize;

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

fn main() {
    let cli = Cli::parse();
    let content = std::fs::read_to_string(cli.config).unwrap();
    let config: Config = toml::from_str(&content).unwrap();

    if let Some(real) = config.real {
        let runs = filter_available_runs(&real);
        write_samplesheet(runs);
    } else {
        println!("No real patients...");
    }

    if let Some(silico) = config.silico {
        println!("Silico patients: {:?}", silico.patients);
    } else {
        println!("No silico patients...");
    }
    match &cli.command {
        Some(Commands::Samplesheet {}) => {}
        None => {}
    }

    // Continued program logic goes here...
}
