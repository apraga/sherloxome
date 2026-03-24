pub mod analyze;
pub mod cli;
pub mod controls;
pub mod fastqbaid2020;
pub mod giab;
pub mod plot;
use log;
use std::path::PathBuf;

/// Helper to download a single URL. Assume the output directory exist
fn download_blocking(url: &str, out: &PathBuf) {
    if out.exists() {
        log::info!("{:?} already exists, skipping", out);
    } else {
        let resp = reqwest::blocking::get(url)
            .expect("Failed to download")
            .bytes()
            .expect("Invalid body in download");
        std::fs::write(out, resp).expect("Failed to write dowloaded file");
        log::debug!("Downloaded");
    }
}
