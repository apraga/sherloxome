use clap::Parser;
use glob::glob;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
// use std::io::{Error, ErroKind};
use std::error::Error;

fn capture(path: &Path) -> Result<PathBuf, Box<dyn Error>> {
    let path_str = path.to_str().ok_or("Path is not valid UTF-8")?;
    let capt_base = PathBuf::from("../../baid2020/bed");

    let capt_file = if path_str.contains("idt") {
        "idt_capture.grch38.bed"
    } else if path_str.contains("truseq") {
        "truseq-dna-exome-targeted-regions-manifest-v1-2-hg38.bed"
    } else {
        println!("Agilent by default");
        "agilent.targets.grch38.bed"
    };

    Ok(capt_base.join(capt_file))
}

fn mosdepth(path: PathBuf, outdir: &Path) -> Result<(), Box<dyn Error>> {
    let stem = path.file_stem().ok_or("Fails to extract file stem")?.to_str().ok_or("Stem is not valid UTF-8")?;
    let prefix = outdir.join(stem);
    let fasta = PathBuf::from("../../genome/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz");
    let capt = capture(&path)?;
    let capt_str = capt.to_str().ok_or("Capture path is not valid UTF-8")?;

    assert!(fasta.exists());
    assert!(capt.exists());

    println!("{:?} : {}", stem, capt.display());
    let output = Command::new("mosdepth")
        .args(["-x", "--threads", "6", "--by", capt_str, "--fasta", fasta.to_str().unwrap()])
        .args([prefix.to_str().unwrap(), path.to_str().unwrap()])
        .output()?;

    println!("{}", String::from_utf8(output.stdout)?);
    println!("{}", String::from_utf8(output.stderr)?);
    Ok(())
}

fn coverage(indir: &str, outdir: &Path) -> Result<(), Box<dyn Error>> {
    let expr = format!("{}/**/HG*.bam", indir);
    for entry in glob(&expr)? {
        let path = entry?;
        mosdepth(path, outdir)?;
    }
    Ok(())
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Directory searched for bam files (recursive)
    #[arg(short, long)]
    dir: String,

    /// Output directory
    #[arg(short, long)]
    out: String,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let outdir = PathBuf::from(&args.out);
    fs::create_dir_all(&outdir)?;

    coverage(&args.dir, &outdir)?;

    Ok(())
}
