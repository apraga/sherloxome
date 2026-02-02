use clap::{Parser, Subcommand};
use sherloxome::giab::Patient;
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

fn main() {
    let cli = Cli::parse();

    println!("Value for config: {}", cli.config.display());

    let content = std::fs::read_to_string(cli.config).unwrap();
    let config: Config = toml::from_str(&content).unwrap();
    println!("{:?}", config);
    if let Some(real) = config.real {
        println!("Real patients: {:?}", real.patients);
    } else {
        println!("No real patients...");
    }

    if let Some(silico) = config.silico {
        println!("Silico patients: {:?}", silico.patients);
    } else {
        println!("No silico patients...");
    }
    // if let Some(config_path) = cli.config.as_deref() {
    // println!("Value for config: {}", config_path.display());
    // }

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match &cli.command {
        Some(Commands::Samplesheet {}) => {}
        None => {}
    }

    // Continued program logic goes here...
}
