use std::process::Command;
use glob::glob;
use clap::Parser;
use std::path::{Path,PathBuf};
use std::fs;
use flate2::read::GzDecoder;
use std::io;

fn query_bed(path: &Path) -> PathBuf {
    let path_str = path.to_str().expect("Path is not valid UTF-8").to_lowercase();
    let capt_base = PathBuf::from("../../baid2020/bed");

    let capt_file = if path_str.contains("idt") {
        "idt_capture.grch38.bed"
    } else if path_str.contains("truseq") {
        "truseq-dna-exome-targeted-regions-manifest-v1-2-hg38.bed"
    } else {
        println!("Agilent by default");
        "agilent.targets.grch38.bed"
    };

     capt_base.join(capt_file)
}


fn truth_base() -> PathBuf {
   PathBuf::from("../../baid2020/giab")
} 

fn truth_vcf(patient: &str) -> PathBuf {
    let vcf = PathBuf::from(format!("{patient}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"));
    truth_base().join(vcf)
}

fn truth_bed(patient: &str) -> PathBuf {
    let bed = PathBuf::from(format!("{patient}_GRCh38_1_22_v4.2.1_benchmark.bed"));
    truth_base().join(bed)
}

fn patient(vcf: &PathBuf) -> String {
    let vcf_ = vcf.to_str().expect("Fails to convert VCF to UTF8");
    if vcf_.contains("HG001") { String::from("HG001")}
    else if vcf_.contains("HG002") { String::from("HG002")}
    else if vcf_.contains("HG003") { String::from("HG003")}
    else if vcf_.contains("HG004") { String::from("HG004")}
    else if vcf_.contains("HG005") { String::from("HG005")}
    else if vcf_.contains("HG006") { String::from("HG006")}
    else if vcf_.contains("HG007") { String::from("HG007")}
    else { panic!("No patient found")}
}

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

fn happy_vcfeval(truth_vcf: &PathBuf, truth_bed: &PathBuf, query_vcf: &PathBuf, query_bed: &PathBuf, 
                    out: &PathBuf, fasta: &PathBuf, sdf: &PathBuf) {
    let summary = out.join(".summary.csv");
    if summary.exists() {
        println!("Summary file already exists (summary)");
    } else {
        let args = ["hap.py", 
                    truth_vcf.to_str().unwrap(), 
                    query_vcf.to_str().unwrap(), 
                    "-o", out.to_str().unwrap(), 
                    "--false-positives", truth_bed.to_str().unwrap(), 
                    "--target-regions", query_bed.to_str().unwrap(), 
                    "--reference", fasta.to_str().unwrap(),
                    "--engine=vcfeval",
                    "--engine-vcfeval-template", sdf.to_str().unwrap()];
        let output = Command::new("hap.py")
            .args(&args)
            .output()
            .expect("Failed to run hap.py from rust");
        if !output.status.success() {
            println!("Happy fails to execute:");
            println!("{:?}", String::from_utf8(output.stdout))
        }
    }
}

fn extract_fasta(fasta_gz: &PathBuf) -> PathBuf {
    let fasta = fasta_gz.with_extension("");
    if !fasta.exists() {
        let gz_file = fs::File::open(fasta_gz)
            .expect("Fails to read gzip fasta");
        let mut gz_decoder = GzDecoder::new(gz_file);
        let mut file = fs::File::create(fasta.clone())
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


fn compare(indir: &str, outdir: &PathBuf) {
    let genome = PathBuf::from("../../genome/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz");
    assert!(genome.exists(), "Missing genome {}", genome.display());
    let fasta = extract_fasta(&genome);
    let sdf = fasta_to_sdf(&fasta);

    let expr = format!("{}/**/HG*.vcf", indir);
    for entry in glob(&expr).expect("Fails to find vcf") {
        compare_vcf(&entry.unwrap(), &outdir, &fasta, &sdf);
    }
}

fn compare_vcf(entry: &PathBuf, outdir: &PathBuf, fasta: &PathBuf, sdf: &PathBuf) {
    let patient = patient(&entry);
    let vcf_truth = truth_vcf(&patient);
    assert!(vcf_truth.exists(), "Missing truth vcf {}", vcf_truth.display());
    let bed_truth = truth_bed(&patient);
    assert!(bed_truth.exists(), "Missing truth bed {}", bed_truth.display());
    let bed_query = query_bed(&entry);
    assert!(bed_query.exists(), "Missing query bed {}", bed_query.display());
    happy_vcfeval(&vcf_truth, &bed_truth, &entry, &bed_query, &outdir, &fasta, &sdf);
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Directory searched for vcf files (recursive)
    #[arg(short, long)]
    dir: String,

    /// Output directory
    #[arg(short, long)]
    out: String,
}


fn main() {
    let args = Args::parse();

    let outdir = PathBuf::from(&args.out);
    fs::create_dir_all(&outdir).expect("Failed to create output dir");

    compare(&args.dir, &outdir);
    // println!("{}", genome.exists());
    // println!("{}", query_vcf("hiseq4000", "agilent", "50x", "HG001", "gatk4").display());
}
