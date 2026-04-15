pub mod analyze;
use std::time::Duration;
pub mod cli;
pub mod fastqbaid2020;
pub mod giab;
pub mod plot;
pub mod setup;
pub mod silico;
use log;
use std::path::PathBuf;

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
