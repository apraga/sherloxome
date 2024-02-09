use std::process::Command;
use std::path::{Path,PathBuf};
use std::fs::File;
use flate2::read::GzDecoder;
use std::io;

fn fasta_to_sdf(genome: &PathBuf) -> PathBuf {
    let sdf = genome.with_extension("sdf");
    if !sdf.exists() {
        Command::new("rtg")
            .arg("format")
            .arg(genome)
            .arg("-o")
            .arg(&sdf)
            .status()
            .expect("Failed to execute command");
    }
    else {
        println!("{} already exists", sdf.display())
    }
    sdf
}

fn happy_vcfeval(truth_vcf: &str, truth_bed: &str, query_vcf: &str, query_bed: &str, outdir: &str, fasta: &str, fai: &str, sdf: &str) {
    let summary = format!("{}.summary.csv", outdir);
    if Path::new(&summary).exists() {
        println!("Summary file already exists (summary)");
    } else {
        let args = ["hap.py", truth_vcf, query_vcf, "-o", outdir, "--false-positives", truth_bed, 
                "--target-regions", query_bed, "--reference", fasta,
                "--engine=vcfeval",
                "--engine-vcfeval-template", sdf];
        let output = Command::new("hap.py")
            .args(&args)
            .output()
            .expect("Failed to execute command");
    }
}

fn query_bed(capture: &str) -> PathBuf {
    match capture  {
        "Agilent v7" => PathBuf::from( "../baid2020/bed/agilent.targets.grch38.bed"),
        "TruSeq" => PathBuf::from("../baid2020/bed/grch38_refseq.bed"),
        "IDT-xGen" => PathBuf::from("../baid2020/bed/idt_capture.grch38.bed"),
        _ => {
            println!("No BED for capture type {} !", capture);
            std::process::exit(1);
        }
    }
}

fn query_vcf(sequencer: &str, capture: &str, coverage: &str, patient: &str, caller: &str) -> PathBuf {
    let root = PathBuf::from("../baid2020/grch38/vcf")
        .join(sequencer)
        .join(format!("wes_{}", capture))
        .join(coverage);
    let vcf =format!("{}.{}.wes-{}.{}.{}.grch38.vcf.gz", patient, sequencer, capture, coverage, caller);
    root.join(vcf)
}

fn extract_fasta(fasta_gz: &PathBuf) -> PathBuf {
    let fasta = fasta_gz.with_extension("");
    if !fasta.exists() {
        let gz_file = File::open(fasta_gz)
            .expect("Fails to read gzip fasta");
        let mut gz_decoder = GzDecoder::new(gz_file);
        let mut file = File::create(fasta.clone())
            .expect("Fails to create fasta");
        io::copy(&mut gz_decoder, &mut file)
            .expect("failed to extract fasta");
        println!("Extracted {}", fasta.display());
    }
    else {
        println!("Already extracted {}", fasta.display());
    }
    fasta
}

fn main() {
    println!("Hi");
    let genome = PathBuf::from("../../genome/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz");
    let fasta = extract_fasta(&genome);
    let sdf = fasta_to_sdf(&fasta);
    // println!("{}", genome.exists());
    // println!("{}", query_vcf("hiseq4000", "agilent", "50x", "HG001", "gatk4").display());
}
