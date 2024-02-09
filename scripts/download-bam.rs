use std::fs;
use std::process::Command;

fn generate_url(sequencer: &str, capture: &str, coverage: &str, patient: &str) -> String {
    let root = "https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam";
    format!("{root}/{sequencer}/wes_{capture}/{coverage}/{patient}.{sequencer}.wes-{capture}.{coverage}.dedup.grch38.bam")
}

fn download(sequencer: &str, capture: &str, coverage: &str) {
    let outdir = format!("../baid2020/grch38/bam/{sequencer}/wes_{capture}/{coverage}");
    fs::create_dir_all(&outdir).expect("Failed to create output directory");

    let patients: Vec<_> = (1..=7).map(|x| format!("HG00{x}")).collect();
    let urls: Vec<_> = patients.iter()
                               .map(|patient| generate_url(sequencer, capture, coverage, patient))
                               .collect();

    for url in urls {
        println!("Downloading {} {} : {}", sequencer, capture, url);
        Command::new("sdf")
                .arg("get")
                .arg(&url)
                .output()
                .expect("Failed to execute sdf command");
    }
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
