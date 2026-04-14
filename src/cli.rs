use crate::analyze::analyze;
use crate::controls::{check_controls_deps, edit_bam, index_bam, remove_hard_clips, resolve_bam};
use crate::download_blocking;
use crate::fastqbaid2020::{Capture, Depth, Run, Sequencer};
use crate::fastqbaid2020::{available, samplesheet_real};
use crate::giab::{Patient, bed_url, vcf_url};
use crate::giab::{all_patients, bed_file, vcf_file};
use crate::plot::plot;
use clap::{Parser, Subcommand};
use itertools::iproduct;
// use polars::io::resolve_homedir;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::path::PathBuf;
use toml;

use serde::Deserialize;

/// Configuration file for the user definig which GIAB data and which in silico data
#[derive(Deserialize, Debug)]
struct Config {
    fasta: PathBuf,
    real: Option<RealConfig>,
    silico: Option<SilicoConfig>,
    capture: HashMap<String, String>,
}

#[derive(Deserialize, Debug)]
struct SilicoConfig {
    capture: Capture,
    // fastq: bool,
    fasta: PathBuf,
    bam: Option<String>,
}

#[derive(Deserialize, Debug)]
struct RealConfig {
    patients: Vec<Patient>,
    sequencers: Vec<Sequencer>,
    captures: Vec<Capture>,
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
        /// Sets a custom config file (TOML) for each command
        #[arg(short, long, value_name = "FILE", default_value = "config.toml")]
        config: PathBuf,
    },
    /// Analyze VCF with hap.py
    Analyze {
        /// Sets a custom config file (TOML) for each command
        #[arg(short, long, value_name = "FILE", default_value = "config.toml")]
        config: PathBuf,

        #[arg(short, long, value_name = "INPUT_DIR")]
        input: PathBuf,

        #[arg(short, long, value_name = "OUTPUT_DIR")]
        output: PathBuf,
    },
    /// Generate controls from clinvar data
    Controls {
        /// Sets a custom config file (TOML) for each command
        #[arg(short, long, value_name = "FILE", default_value = "config.toml")]
        config: PathBuf,
    },
    Download {
        #[arg(short, long, value_name = "OUTPUT_DIR", default_value = "data/ref")]
        output: PathBuf,
    },
    Plot {
        #[arg(short, long, value_name = "INPUT_FILE")]
        input: PathBuf,
        #[arg(short, long, value_name = "OUTPUT_FILE", default_value = "plot.html")]
        output: PathBuf,
    },
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

/// Generate a samplesheet for real and in-silico data according to a configuration file
fn generate_samplesheet(config_file: PathBuf) {
    let config = read_config(config_file);
    if let Some(real) = &config.real {
        let runs = filter_available_runs(&real);
        samplesheet_real(runs);
    } else {
        println!("No real patients...");
    }

    if let Some(silico) = &config.silico {
        samplesheet_silico(silico);
    }
}

fn samplesheet_silico(config: &SilicoConfig) {}

/// Download all reference VCF and BED in a directory.
fn download(out_dir: PathBuf) {
    std::fs::create_dir_all(&out_dir).expect("Could not create directory");
    for p in all_patients() {
        println!("{:?}", p);
        let mut vcf = out_dir.clone().join(vcf_file(&p));
        let bed = out_dir.clone().join(bed_file(&p));
        download_blocking(&vcf_url(&p), &vcf);
        vcf.set_extension("gz.tbi");
        download_blocking(&vcf_url(&p), &vcf);
        download_blocking(&bed_url(&p), &bed);
    }
}

/// Wrapper to generate controls from clinvar data for several captures files
fn generate_controls(conf: Config) {
    check_controls_deps();
    if let Err(e) = try_generate_controls(conf) {
        log::error!("Failed to generate insilico controls: {e}");
    }
}

fn try_generate_controls(conf: Config) -> Result<(), Box<dyn Error>> {
    let silico = conf.silico.ok_or("missing [silico] in config")?;
    let capture = silico.capture;
    let outdir = PathBuf::from("data/exp_raw");
    std::fs::create_dir_all(&outdir)?;

    let bed = conf
        .capture
        .get(&capture.to_string())
        .ok_or_else(|| format!("No bed file for {capture}"))?;
    if let Some(bam) = silico.bam {
        log::info!("Editing BAM to insert control : {:?}", bam);
        let bam_path = resolve_bam(bam, &outdir)?;
        index_bam(&bam_path)?;
        let bam_cleaned = remove_hard_clips(&bam_path)?;
        edit_bam(bam_cleaned, bed, capture.to_string(), conf.fasta)?;
    } else {
        log::info!("No BAM editing...");
    }

    Ok(())
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
        Some(Commands::Samplesheet { config }) => {
            log::info!("Generating samplesheet...");
            generate_samplesheet(config.to_path_buf());
        }
        Some(Commands::Analyze {
            config,
            input,
            output,
        }) => {
            log::info!("Analyzing runs...");
            let conf = read_config(config.to_path_buf());
            if let Err(e) = analyze(
                &conf.fasta,
                conf.capture,
                input.to_path_buf(),
                output.to_path_buf(),
            ) {
                eprintln!("Analysis failed: {e}");
                std::process::exit(1);
            }
        }

        Some(Commands::Controls { config }) => {
            log::info!("Generating controls...");
            let conf = read_config(config.to_path_buf());
            generate_controls(conf);
        }

        Some(Commands::Download { output }) => {
            download(output.to_path_buf());
        }
        Some(Commands::Plot { input, output }) => {
            let _ = plot(input.to_path_buf(), output.to_path_buf());
        }

        None => {}
    }
}
