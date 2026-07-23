//! Exome variant-calling benchmarking toolkit.
//!
//! 1. Download GIAB reference VCFs/BEDs and prepare FASTQ data for either real-patient or in silico control, then generates a samplesheet for sarek ([`setup`])
//! 2. Benchmark pipeline output VCFs against GIAB truth sets with hap.py ([`benchmark`])
//! 3. Visualise performances metrics ([`plot`])
pub mod benchmark;
pub mod run;
use flate2::read::GzDecoder;
use std::time::Duration;
use url::Url;
pub mod cli;
use std::error::Error;
pub mod baid2020;
pub mod giab;
pub mod plot;
pub mod setup;
pub mod silico;
pub mod simuscop;
pub mod varben;
use log;
use std::path::PathBuf;
use which::which;

/// Helper to download a single URL. Assume the output directory exist
fn download_blocking(url: &str, out: &PathBuf) {
    if out.exists() {
        log::info!("{:?} already exists, skipping", out);
    } else {
        let client = reqwest::blocking::Client::builder()
            .timeout(None)
            .connect_timeout(Duration::from_secs(30))
            .build()
            .expect("Failed to build HTTP client");

        let mut resp = client
            .get(url)
            .send()
            .expect(&format!("Failed to download {url}"))
            .error_for_status()
            .expect(&format!("HTTP error for {url}"));

        let mut file = std::fs::File::create(out).expect(&format!("Failed to create {out:?}"));

        resp.copy_to(&mut file)
            .expect(&format!("Failed to write to {out:?}"));
        log::debug!("Downloaded to {:?}", out);
    }
}

/// Check if fasta file exists. If not, download GRCh38 full analysis set from NCBI and unzip it.
pub fn resolve_fasta(fasta: &str) -> Result<PathBuf, Box<dyn Error>> {
    let fasta_path = PathBuf::from(fasta);
    let outdir = PathBuf::from("data/ref");
    if fasta.starts_with("http://") || fasta.starts_with("https://") {
        download_file(fasta, &outdir)
    } else if fasta_path.exists() {
        Ok(fasta_path)
    } else {
        log::debug!(
            "{:?} not found, downloading GRCh38 reference from NCBI...",
            fasta
        );
        download_fasta(outdir)
    }
}

fn download_fasta(outdir: PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    let fasta_gz_dist = "GCA_000001405.15_GRCh38_full_analysis_set.fna.gz";
    let fasta_dist = "GCA_000001405.15_GRCh38_full_analysis_set.fna";
    let fai_dist = "GCA_000001405.15_GRCh38_full_analysis_set.fna.fai";
    let root_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines";

    let gz_path = outdir.join(fasta_gz_dist);
    let fai_path = outdir.join(fai_dist);
    let fasta_path = outdir.join(fasta_dist);

    if !fasta_path.exists() {
        let url = format!("{root_url}/{fasta_gz_dist}");
        download_blocking(&url, &gz_path);
        log::debug!("Decompressing {:?}...", gz_path);
        let gz_file = std::fs::File::open(&gz_path).expect("Failed to open gz file");
        let mut decoder = GzDecoder::new(gz_file);
        let mut out_file =
            std::fs::File::create(&fasta_path).expect("Failed to create decompressed fasta");
        std::io::copy(&mut decoder, &mut out_file).expect("Failed to decompress fasta");
    }
    let url = format!("{root_url}/{fai_dist}");
    download_blocking(&url, &fai_path);
    Ok(fasta_path)
}

/// Download a bam file if it's an url, otherwise check the file exist
pub fn resolve_bam(bam: &str, outdir: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    if bam.starts_with("http://") || bam.starts_with("https://") {
        download_file(bam, outdir)
    } else {
        let path = PathBuf::from(bam);
        if !path.exists() {
            return Err(format!("BAM file not found: {}", path.display()).into());
        }
        Ok(path)
    }
}

fn download_file(url: &str, outdir: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    let parsed = Url::parse(&url)?;
    log::info!("Downloading file {}", url);
    let filename = parsed
        .path_segments()
        .and_then(|s| s.last())
        .filter(|s| !s.is_empty())
        .ok_or("could not extract filename from URL")?;

    let output = outdir.join(filename);
    download_blocking(&url, &output);
    Ok(output)
}

/// Panic if dependcies are missing from PATH
pub fn check_deps(tools: &[&str]) {
    log::debug!("Checking dependencies");
    for tool in tools {
        match which(tool) {
            Ok(path) => log::debug!("{tool} ✓  ({})", path.display()),
            Err(_) => {
                log::error!("{tool} ✗  not found in PATH");
                log::error!("Please install dependencies (nix is the recommended way)");
                panic!();
            }
        }
    }
}

fn ref_dir() -> PathBuf {
    PathBuf::from("data/ref")
}
