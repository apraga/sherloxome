//! Command-line interface.
//!
//! Three subcommands:
//! - `setup`   — download GIAB data and generate in silico controls. It is made of steps
//!   (`prepare`, `simuscop`, `varben`, `samplesheet`) that can also be run one by one, for
//!   instance one `simuscop` per profile in a Slurm job array
//! - `benchmark` — benchmark VCF files by comparing to reference VCF
//! - `plot`    — display F1-score boxplots from `merged.csv`
use crate::benchmark::analyze;
use crate::plot::plot;
use crate::setup::{Config, prepare, samplesheet, setup, simuscop, varben};
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
    ///
    /// Without a step, runs `prepare`, every simuscop and varben run of the configuration,
    /// then `samplesheet`.
    Setup {
        /// Sets a custom config file (TOML) for each command
        #[arg(
            short,
            long,
            value_name = "FILE",
            default_value = "config.toml",
            global = true
        )]
        config: PathBuf,

        #[command(subcommand)]
        step: Option<SetupStep>,
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

/// A single step of `setup`, to run them separately (e.g. as a Slurm job array)
#[derive(Subcommand)]
enum SetupStep {
    /// Run once, before the other steps: download reference data and sample clinvar/dbSNP
    ///
    /// Sampling is random and its result is shared by all runs, so it must not happen
    /// in the concurrent `simuscop` and `varben` runs, which fail if it has not been done.
    Prepare {
        /// Simuscop profiles that will be run. Default: the one of the configuration
        #[arg(long, value_name = "PROFILE", num_args = 1..)]
        profile: Vec<PathBuf>,
        /// Varben BAMs that will be run. Default: those of the configuration
        #[arg(long, value_name = "BAM", num_args = 1..)]
        bam: Vec<PathBuf>,
    },
    /// Generate the FASTQ of a single simuscop profile (capture and depth from its filename)
    Simuscop {
        #[arg(long, value_name = "PROFILE")]
        profile: PathBuf,
    },
    /// Insert the sampled variants into a single BAM and generate its FASTQ
    Varben {
        #[arg(long, value_name = "BAM")]
        bam: PathBuf,
    },
    /// Write one samplesheet per capture kit from the FASTQ found in the silico directory
    Samplesheet,
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
        Commands::Setup { config, step } => {
            let conf = read_config(config)?;
            match step {
                None => {
                    log::info!("Setting up runs...");
                    setup(conf)?;
                }
                Some(SetupStep::Prepare { profile, bam }) => prepare(&conf, profile, bam)?,
                Some(SetupStep::Simuscop { profile }) => simuscop(&conf, profile)?,
                Some(SetupStep::Varben { bam }) => varben(&conf, bam)?,
                Some(SetupStep::Samplesheet) => samplesheet(&conf)?,
            }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Cli {
        Cli::try_parse_from(std::iter::once("sherloxome").chain(args.iter().copied())).unwrap()
    }

    #[test]
    fn setup_without_step_runs_everything() {
        let Commands::Setup { step, config } = parse(&["setup"]).command else {
            panic!("expected setup");
        };
        assert!(step.is_none());
        assert_eq!(config, PathBuf::from("config.toml"));
    }

    #[test]
    fn config_is_accepted_after_the_step() {
        let Commands::Setup { step, config } = parse(&[
            "setup",
            "simuscop",
            "--profile",
            "a.profile",
            "-c",
            "x.toml",
        ])
        .command
        else {
            panic!("expected setup");
        };
        assert_eq!(config, PathBuf::from("x.toml"));
        assert!(
            matches!(step, Some(SetupStep::Simuscop { profile }) if profile == PathBuf::from("a.profile"))
        );
    }

    #[test]
    fn prepare_takes_several_profiles_and_bams() {
        let Commands::Setup { step, .. } = parse(&[
            "setup",
            "prepare",
            "--profile",
            "a.profile",
            "b.profile",
            "--bam",
            "c.bam",
        ])
        .command
        else {
            panic!("expected setup");
        };
        let Some(SetupStep::Prepare { profile, bam }) = step else {
            panic!("expected prepare");
        };
        assert_eq!(profile.len(), 2);
        assert_eq!(bam, vec![PathBuf::from("c.bam")]);
    }
}
