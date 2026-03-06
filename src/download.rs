use std::path::Path;
use std::process::Command;

const BASE_URL: &str = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release";

// (sample_id, url_tail, bed_prefbix)
// Exemple of base url  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/
const SAMPLES: &[(&str, &str)] = &[
    ("HG001", "NA12878_HG001", ".bed"),
    (
        "HG002",
        "AshkenazimTrio/HG002_NA24385_son",
        "_noinconsistent.bed",
    ),
    (
        "HG003",
        "AshkenazimTrio/HG003_NA24149_father",
        "_noinconsistent.bed",
    ),
    (
        "HG004",
        "AshkenazimTrio/HG004_NA24143_mother",
        "_noinconsistent.bed",
    ),
    ("HG005", "ChineseTrio/HG005_NA24631_son", ".bed"),
    ("HG006", "ChineseTrio/HG006_NA24694_father", ".bed"),
    ("HG007", "ChineseTrio/HG007_NA24695_mother", ".bed"),
];

// Suffixes shared by every sample; prefixed at runtime with "{id}_"
const FILE_SUFFIXES: &[&str] = &[
    "GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    "GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
];

// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct DownloadResult {
    pub sample: String,
    pub file: String,
    pub status: DownloadStatus,
}

#[derive(Debug)]
pub enum DownloadStatus {
    Downloaded,
    Skipped,
    Failed(String),
}

/// Download all  GIAB truth VCF/BED files into `ref_dir`.
/// Pass an empty `samples` slice to download all seven samples.
pub fn download(ref_dir: &Path) -> anyhow::Result<Vec<DownloadResult>> {
    std::fs::create_dir_all(ref_dir)?;

    let mut results = Vec::new();

    for (id, url_tail) in SAMPLES {
        println!("--- {id} ---");
        for filename in files_for(id) {
            let dest = ref_dir.join(&filename);

            if dest.exists() {
                println!("  [skip] {filename} (already exists)");
                continue;
            }

            let url = format!("{BASE_URL}/{url_tail}/{filename}");
            println!("  [wget] {url}");

            // let result_status = match Command::new("wget")
            //     .args(["--continue", "-P"])
            //     .arg(ref_dir)
            //     .arg(&url)
            //     .status()
            // {
            //     Ok(s) if s.success() => DownloadStatus::Downloaded,
            //     Ok(s) => {
            //         let msg = format!("wget exited with {s}");
            //         eprintln!("  [fail] {filename}: {msg}");
            //         DownloadStatus::Failed(msg)
            //     }
            //     Err(e) => {
            //         let msg = format!("failed to spawn wget: {e}");
            //         eprintln!("  [fail] {filename}: {msg}");
            //         DownloadStatus::Failed(msg)
            //     }
            // };

            // results.push(DownloadResult {
            //     sample: id.to_string(),
            //     file: filename,
            //     status: result_status,
            // });
        }
    }

    let downloaded = results
        .iter()
        .filter(|r| matches!(r.status, DownloadStatus::Downloaded))
        .count();
    let skipped = results
        .iter()
        .filter(|r| matches!(r.status, DownloadStatus::Skipped))
        .count();
    let failed = results
        .iter()
        .filter(|r| matches!(r.status, DownloadStatus::Failed(_)))
        .count();
    println!("\nSummary: {downloaded} downloaded, {skipped} skipped, {failed} failed");

    if failed > 0 {
        anyhow::bail!("{failed} file(s) failed to download");
    }
    Ok(results)
}

/// Print all known sample IDs (useful for shell completion).
pub fn list_samples() {
    for (id, _) in SAMPLES {
        println!("{id}");
    }
}
