use std::process::Command;
use glob::glob;
use std::path::{Path, PathBuf};
use std::str;

fn capture(stem: &str) -> PathBuf {
    if stem.contains("idt") {
        return PathBuf::from("idt_capture.grch38.bed");
    }
    else if stem.contains("truseq") {
        return PathBuf::from("truseq-dna-exome-targeted-regions-manifest-v1-2-hg38.bed");
    }
    else {
        println!("Agilent by default");
        return PathBuf::from("agilent.targets.grch38.bed");
    }

}
fn mosdepth(path: PathBuf) {
    let stem = path.file_stem().expect("Fails to extract file stem").to_str().unwrap();
    // TODO mkdir this directory
    let out = String::from("baid2020/") + stem;
    let dir = path.parent().unwrap();
    let fasta = "../../genome/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz";
    let capt_path = PathBuf::from("../../baid2020/bed").join(capture(stem));
    let capt = capt_path.to_str().unwrap() ;
    println!("{} : {}", stem, capt);
    let output = Command::new("mosdepth")
        .args(["-x", "--threads", "6", "--by",  &capt, "--fasta", fasta])
        .args([&out, path.to_str().unwrap()])
        .output
        ().expect("Fails to run command");
    println!("{}", str::from_utf8(&output.stdout).unwrap());
    println!("{}", str::from_utf8(&output.stderr).unwrap());
}

fn main() {
    for entry in glob("../../baid2020/grch38/bam/**/*.bam").expect("Failed to read glob pattern") {
        match entry {
            Ok(path) => mosdepth(path),
            Err(e) => println!("{:?}", e),
        }
    }
}
