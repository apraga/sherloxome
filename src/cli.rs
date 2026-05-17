//! Command-line interface.
//!
//! Three subcommands:
//! - `setup`   — download GIAB data and generate in silico controls
//! - `benchmark` — benchmark VCF files by comparing to reference VCF
//! - `plot`    — display F1-score boxplots from `merged.csv`
use crate::benchmark::analyze;
use crate::plot::plot;
use crate::setup::{Config, setup};
use clap::{Parser, Subcommand};
use std::path::PathBuf;
use toml;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Download or generate FASTQ and generate samplesheet
    Setup {
        /// Sets a custom config file (TOML) for each command
        #[arg(short, long, value_name = "FILE", default_value = "config.toml")]
        config: PathBuf,
    },
    /// Benchmark all VCF in a directory with hap.py
    Benchmark {
        /// Sets a custom config file (TOML) for each command
        #[arg(short, long, value_name = "FILE", default_value = "config.toml")]
        config: PathBuf,

        #[arg(short, long, value_name = "INPUT_DIR")]
        input: PathBuf,

        #[arg(short, long, value_name = "OUTPUT_DIR")]
        output: PathBuf,
    },
    Plot {
        #[arg(short, long, value_name = "INPUT_FILE")]
        input: PathBuf,
        #[arg(short, long, value_name = "OUTPUT_FILE", default_value = "plot.html")]
        output: PathBuf,
    },
}

fn read_config(fname: PathBuf) -> Config {
    let content = std::fs::read_to_string(fname).unwrap();
    let conf: Config = toml::from_str(&content).unwrap();
    conf
}

/// Read CLI arguments and call relevant functions
pub fn process_cli() {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Setup { config } => {
            let conf = read_config(config.to_path_buf());
            log::info!("Setting up runs...");
            if let Err(e) = setup(conf) {
                eprintln!("Setup failed: {e}");
                std::process::exit(1);
            }
        }
        Commands::Benchmark {
            config,
            input,
            output,
        } => {
            let conf = read_config(config.to_path_buf());
            log::info!("Analyzing runs...");
            if let Err(e) = analyze(&conf, input.to_path_buf(), output.to_path_buf()) {
                eprintln!("Analysis failed: {e}");
                std::process::exit(1);
            }
        }
        Commands::Plot { input, output } => {
            let _ = plot(input.to_path_buf(), output.to_path_buf());
        }
    }
}
