use std::fs;
use std::process::Command;
use std::io::{self, Write};
use std::path::PathBuf;



// Return the url down to the directory level for a sets of patients
fn root_url(sequencer: &str, capture: &str, coverage: &str) -> String {
    let root = "https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam";
    format!("{root}/{sequencer}/wes_{capture}/{coverage}")
}

// Retur nthe bam filename
fn bam_file(sequencer: &str, capture: &str, coverage: &str, patient: &str) -> String {
    // Underscore in filename for idt
    if capture == "idt" || 
        (sequencer == "novaseq" && capture == "truseq" && patient != "HG006"){ 
        format!("{patient}.{sequencer}.wes_{capture}.{coverage}.dedup.bam.bai") 
    }
    else { 
        format!("{patient}.{sequencer}.wes-{capture}.{coverage}.dedup.grch38.bam.bai") 
    }
}

fn download(sequencer: &str, capture: &str, coverage: &str) {
    let outdir = format!("../baid2020/grch38/bam/{sequencer}/wes_{capture}/{coverage}");
    fs::create_dir_all(&outdir).expect("Failed to create output directory");

    let patients: Vec<_> = (1..=7).map(|x| format!("HG00{x}")).collect();
    for patient in patients {
        let dir = root_url(&sequencer, &capture, &coverage);
        let bam = bam_file(&sequencer, &capture, &coverage, &patient);
        let url = format!("{dir}/{bam}");
        if PathBuf::from(&outdir).join(&bam).exists() {
            println!("{bam} already exists, skipping");
        }
        else {
            download_url(&url, &outdir)
        }
    }
}

fn download_url(url: &str, outdir: &str) {
    println!("Downloading {url} into {outdir}");
    let output = Command::new("sdf")
                    .arg("get")
                    .arg(&url)
                    .current_dir(&outdir)
                    .output()
                    .expect("Failed to execute sdf command");
    io::stderr().write_all(&output.stderr).unwrap();
}

fn main() {
    for sequencer in vec!["hiseq4000", "novaseq"] {
        for capture in  vec!["agilent", "idt", "truseq"] {
            for coverage in vec!["50x"] {
                download(sequencer, capture, coverage);
            }
        }
    }
}
