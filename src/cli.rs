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
use std::error::Error;
use std::path::PathBuf;

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

fn read_config(fname: &PathBuf) -> Result<Config, Box<dyn Error>> {
    let content = std::fs::read_to_string(fname)
        .map_err(|e| format!("Cannot read config {}: {e}", fname.display()))?;
    let conf: Config =
        toml::from_str(&content).map_err(|e| format!("Invalid config {}: {e}", fname.display()))?;
    conf.validate()?;
    Ok(conf)
}

/// Read CLI arguments and call subfunctions
pub fn process_cli() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Setup { config } => {
            let conf = read_config(config)?;
            log::info!("Setting up runs...");
            setup(conf)?;
        }
        Commands::Benchmark {
            config,
            input,
            output,
        } => {
            let conf = read_config(config)?;
            log::info!("Analyzing runs...");
            analyze(&conf, input.clone(), output.clone())?;
        }
        Commands::Plot { input, output } => {
            plot(input.clone(), output.clone())?;
        }
    }
    Ok(())
}
