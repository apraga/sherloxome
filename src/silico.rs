//! # Insilico Controls
//! Generate control variants by sampling clinvar data
use crate::download_blocking;
use flate2::read::GzDecoder;
use log;
use noodles::bed;
use noodles::bgzf;
use noodles::fasta;
use noodles::vcf;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::io::Write as VCFWrite;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::info::field::Value;
use noodles::vcf::variant::record_buf::AlternateBases as AltBasesBuf;
use rand::RngExt;
use rand::prelude::IteratorRandom;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitStatus;
use std::process::{Command, Stdio};
use std::thread;
use url::Url;
use which::which;

pub fn check_controls_deps() {
    let tools = ["bwa", "samtools", "muteditor"];

    log::debug!("Checking dependencies");
    for tool in tools {
        match which(tool) {
            Ok(path) => log::debug!("{tool} ✓  ({})", path.display()),
            Err(_) => panic!("{tool} ✗  not found in PATH"),
        }
    }
}
fn download_clinvar() -> PathBuf {
    log::debug!("Downloading clinvar vcf...");
    let mut url: String =
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz".to_string();
    let vcf = PathBuf::from("data/exp_raw/clinvar.vcf.gz");
    download_blocking(&url, &vcf);

    url.push_str(".tbi");
    let tbi = PathBuf::from("data/exp_raw/clinvar.vcf.gz.tbi");
    download_blocking(&url, &tbi);
    vcf
}

/// Pseudo-interval tree..
type BedIndex = HashMap<String, Vec<(usize, usize)>>;

/// For each chromosome, store sorted list of (start, end) 0-based half-open intervals
/// Assubme bed is stored
fn load_bed(bed: &PathBuf) -> Result<BedIndex, Box<dyn Error>> {
    let mut reader = File::open(bed)
        .map(BufReader::new)
        .map(bed::io::Reader::<3, _>::new)?;

    let mut index = BedIndex::new();
    let mut record = bed::Record::default();

    while reader.read_record(&mut record)? != 0 {
        let chrom_raw = std::str::from_utf8(record.reference_sequence_name())?;
        let chrom = chrom_raw
            .strip_prefix("chr")
            .unwrap_or(chrom_raw)
            .to_string();
        let start = usize::from(record.feature_start()?);
        let end = usize::from(record.feature_end().ok_or("missing BED end position")??);
        index.entry(chrom).or_default().push((start, end));
    }
    Ok(index)
}

/// Is a 1-based VCF position inside any BED interval on this chromosome?
fn in_capture(bed: &BedIndex, chrom: &str, pos: u64) -> bool {
    // Convert 1-based VCF pos to 0-based
    let pos0 = pos as usize - 1;

    let Some(intervals) = bed.get(chrom) else {
        return false;
    };

    // Binary search: find the last interval that starts at or before pos0
    let i = intervals.partition_point(|(start, _)| *start <= pos0);

    // Check if pos0 falls before the end of that interval
    i.checked_sub(1)
        .map(|j| intervals[j].1 > pos0)
        .unwrap_or(false)
}

/// Pathogenic variant are preferred
fn clnsig_priority(sig: &str) -> Option<u8> {
    match sig {
        "Pathogenic" => Some(0),
        "Likely_pathogenic" => Some(1),
        "Uncertain_significance" => Some(2),
        _ => None,
    }
}
/// Only keep VOUS, Pathogenic or Likely Pathogenic or VOUS variant (whatever the number of submisson)
/// If there are several variants, only keep the most severe
/// Info field is complicated... Some(Ok(Some(Array([Ok(Some("Uncertain_significance"))]))))
fn is_not_benign(info: &vcf::record::Info, header: &vcf::Header) -> bool {
    info.get(header, "CLNSIG")
        .and_then(|r| r.ok())
        .flatten()
        .is_some_and(|v| match v {
            Value::String(s) => clnsig_priority(s.as_ref()).is_some(),
            Value::Array(arr) => match arr {
                noodles::vcf::variant::record::info::field::value::Array::String(arr) => arr
                    .iter()
                    .filter_map(|s| s.ok().flatten())
                    .any(|s| clnsig_priority(s.as_ref()).is_some()),
                _ => false,
            },
            _ => false,
        })
}

/// Only keep SNV
fn is_snv(record: &vcf::Record) -> bool {
    record.reference_bases().len() == 1
        && record
            .alternate_bases()
            .iter()
            .filter_map(|a| a.ok())
            .all(|a| a.len() == 1)
}

/// Write VCF record to varben format. Assume chromosoes do NOT have a chr prefix
/// Generate AF betwen 0.4 and 0.6
fn write_varben(
    w: &mut impl Write,
    record: &RecordBuf,
    rng: &mut impl RngExt,
) -> Result<(), Box<dyn Error>> {
    let af: f64 = rng.random_range(0.4..0.6);
    let chrom_raw = record.reference_sequence_name();
    let chrom_nb = chrom_raw.strip_prefix("chr").unwrap_or(chrom_raw);
    let chrom_chr = format!("chr{}", chrom_nb);
    let pos = record.variant_start().unwrap();
    let alt = record.alternate_bases();
    let alt_merged = alt.iter().filter_map(|a| a.ok()).next().unwrap_or(".");
    writeln!(
        w,
        "{}\t{}\t{}\t{}\tsnv\t{}",
        chrom_chr, pos, pos, af, alt_merged
    )?;
    Ok(())
}

/// Sort chromsome by natural ornder
fn chrom_order(chrom: &str) -> (u8, u32) {
    match chrom.trim_start_matches("chr") {
        "X" => (1, 0),
        "Y" => (2, 0),
        "M" | "MT" => (3, 0),
        n => (0, n.parse().unwrap_or(u32::MAX)),
    }
}

/// Sample clinvar variant and prefix "chr" to chromosome names
fn sample_clinvar_variants(
    reader: &mut vcf::io::Reader<bgzf::io::Reader<File>>,
    header: &vcf::Header,
    capture: &BedIndex,
    spacing: u32,
    n: u32,
    rng: &mut impl RngExt,
) -> Vec<RecordBuf> {
    let mut last_pos: HashMap<String, u64> = HashMap::new();
    let mut selected: Vec<RecordBuf> = reader
        .records()
        .filter_map(|r| r.ok())
        // Filter first as we need INFO field
        .filter(|record| keep_variant(record, header, capture, &mut last_pos, spacing))
        // The add chr prefilx
        .map(|record| add_chr_prefix(&record))
        .sample(
            rng,
            n.try_into().expect("Fails to convert nb variants to usize"),
        );
    sort_by_chromosome(&mut selected);
    log::debug!("Selected {} variants", selected.len());
    selected
}

fn sort_by_chromosome(variants: &mut Vec<RecordBuf>) {
    variants.sort_by(|a, b| {
        chrom_order(a.reference_sequence_name())
            .cmp(&chrom_order(b.reference_sequence_name()))
            .then_with(|| {
                let ap = a.variant_start().map(usize::from).unwrap_or(0);
                let bp = b.variant_start().map(usize::from).unwrap_or(0);
                ap.cmp(&bp)
            })
    });
}

/// Select n clinvar pathogenic SNV inside the capture kit 50bp apart.
/// Output is both a mutation file for varben and a VCF (hardcoded into data/exp_raw/clinvar_$CAPTURE)
/// The most efficient way is to parse clinvar once to get variant in the capture kit and 50bp apart.
/// In a second pass, sample randomly n of them.
pub fn sample_clinvar(
    clinvar_vcf: PathBuf,
    bed: PathBuf,
    spacing: u32,
    n: u32,
    vcf_out: PathBuf,
    mut_out: PathBuf,
) -> Result<vcf::Header, Box<dyn Error>> {
    let capture = load_bed(&bed)?;
    let mut reader = File::open(&clinvar_vcf)
        .map(bgzf::io::Reader::new)
        .map(vcf::io::Reader::new)
        .expect("Failed to open clinvar file");
    let header = reader.read_header()?;

    if vcf_out.exists() && mut_out.exists() {
        log::debug!(
            "Skip sampling clinvar as output files already exists: {:?} and {:?}",
            vcf_out,
            mut_out
        );
    } else {
        log::debug!("Sampling {n} variants for insertion");
        let mut rng = rand::rng();
        let variants =
            sample_clinvar_variants(&mut reader, &header, &capture, spacing, n, &mut rng);

        write_sampled_clinvar_vcf(&variants, &header, &vcf_out)?;
        write_sampled_clinvar_mut(&variants, &mut_out, &mut rng)?;
        log::debug!("Wrote {:?} and {:?}", vcf_out, mut_out);
    }
    Ok(header)
}

/// Assume there is no chr prefix for chromosome
fn write_sampled_clinvar_vcf(
    variants: &Vec<RecordBuf>,
    header: &vcf::Header,
    vcf_out: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let mut writer = File::create(&vcf_out)
        .map(bgzf::io::Writer::new)
        .map(vcf::io::Writer::new)?;
    writer.write_header(&header)?;

    for record in variants {
        writer.write_variant_record(header, record)?;
    }
    Ok(())
}

/// Write varben output as VCF successful insertions only, failures are not properly formatted,
pub fn write_varben_as_vcf(
    bam: &PathBuf,
    header: vcf::Header,
    outdir: PathBuf,
    fasta: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let varben_dir = outdir.join("varben");
    let bam_stem = bam.file_stem().unwrap().to_string_lossy();

    let success_list = varben_dir.join("success_list.txt");
    let success_vcf = outdir.join(format!("{bam_stem}_success.vcf.gz"));
    write_varben_as_vcf_single(success_list, success_vcf, header, fasta)
}

/// Helper to write a .vcf.gz
// fn vcf_writer(
//     path: &PathBuf,
//     header: &vcf::Header,
// ) -> Result<vcf::io::Writer<bgzf::io::Writer<File>>, Box<dyn Error>> {
//     let mut writer = File::create(path)
//         .map(bgzf::io::Writer::new)
//         .map(vcf::io::Writer::new)?;
//     writer.write_header(header)?;
//     Ok(writer)
// }

fn write_sampled_clinvar_mut(
    variants: &Vec<RecordBuf>,
    mut_out: &PathBuf,
    rng: &mut impl RngExt,
) -> Result<(), Box<dyn Error>> {
    // Prepare mutation file
    let mut varben_w = BufWriter::new(File::create(&mut_out)?);
    writeln!(varben_w, "#chrom\tstart\tend\tAF\ttype\talt")?;

    for v in variants {
        write_varben(&mut varben_w, v, rng)?;
    }
    varben_w.flush()?;
    Ok(())
}

/// We cannot edit a record directly, recreate it without INFO and sample.
/// This is just for clinvar data for quick check
fn add_chr_prefix(record: &vcf::Record) -> RecordBuf {
    let chrom_raw = record.reference_sequence_name();
    let chrom_nb = chrom_raw.strip_prefix("chr").unwrap_or(chrom_raw);

    let alts: Vec<String> = record
        .alternate_bases()
        .iter()
        .filter_map(|a| a.ok())
        .map(|a| a.to_string())
        .collect();

    let mut builder = RecordBuf::builder()
        .set_reference_sequence_name(format!("chr{}", chrom_nb))
        .set_reference_bases(record.reference_bases().to_string())
        .set_alternate_bases(AltBasesBuf::from(alts));

    if let Some(Ok(pos)) = record.variant_start() {
        builder = builder.set_variant_start(pos);
    }

    builder.build()
}

/// For varian in the capture file and SNV and not bening, save it for writing
/// Update lats saved position (`last_pos)`
/// Input VCF may use chr prefix
fn keep_variant(
    record: &vcf::Record,
    header: &vcf::Header,
    capture: &BedIndex,
    last_pos: &mut HashMap<String, u64>,
    spacing: u32,
) -> bool {
    let Some(Ok(start)) = record.variant_start() else {
        return false;
    };
    let pos = usize::from(start) as u64;
    let chrom_raw = record.reference_sequence_name();
    let chrom = chrom_raw.strip_prefix("chr").unwrap_or(chrom_raw);
    let last = last_pos.get(chrom).copied().unwrap_or(0);

    let ok = in_capture(capture, chrom, pos)
        && is_not_benign(&record.info(), header)
        && is_snv(record)
        && pos != last
        && pos.abs_diff(last) >= spacing.into();

    if ok {
        last_pos.insert(chrom.to_string(), pos);
    }
    ok
}

/// Helper function to index a BAM file
pub fn index_bam(bam: &PathBuf) -> Result<(), Box<dyn Error>> {
    let bai = bam.with_extension("bam.bai");
    // Index the BAM
    if !bai.exists() {
        log::debug!("Indexing bam");
        let status = Command::new("samtools")
            .args(["index", &bam.to_string_lossy()])
            .status()
            .expect("Failed to run samtools index");
        if !status.success() {
            return Err(format!("samtools index exited with status {status}").into());
        }
    } else {
        log::debug!("{:?} already exists ", bai);
    }
    Ok(())
}

/// Bam to fastq conversion require hard clips removal
pub fn remove_hard_clips(bam: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    let stem = bam.file_stem().unwrap().to_string_lossy();
    let parent = bam.parent().unwrap_or(Path::new("."));
    let filtered = parent.join(format!("{stem}_nohardclip.bam"));

    if !filtered.exists() {
        log::info!("Removing hard clips");
        let cmd = format!(
            "samtools view -h {} | awk '$6 !~ /H/{{print}}' | samtools view -bS - > {}",
            bam.display(),
            filtered.display()
        );
        let status = Command::new("sh")
            .args(["-c", &cmd])
            .status()
            .expect("Failed to run hard clip filter");

        if !status.success() {
            return Err(format!("samtools exited with status {status}").into());
        }
    } else {
        log::info!("Hardlink removal already done");
    }
    index_bam(&filtered)?;
    Ok(filtered)
}

/// Check if fasta file exists. If not, download GRCh38 full analysis set from NCBI and unzip it.
pub fn resolve_fasta(fasta: &PathBuf) -> PathBuf {
    if fasta.exists() {
        return fasta.clone();
    }
    log::debug!(
        "{:?} not found, downloading GRCh38 reference from NCBI...",
        fasta
    );
    let fasta_gz_dist = "GCA_000001405.15_GRCh38_full_analysis_set.fna.gz";
    let fasta_dist = "GCA_000001405.15_GRCh38_full_analysis_set.fna";
    let fai_dist = "GCA_000001405.15_GRCh38_full_analysis_set.fna.fai";
    let root_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines";

    let gz_path = PathBuf::from("data/ref").join(fasta_gz_dist);
    let fai_path = PathBuf::from("data/ref").join(fai_dist);
    let fasta_path = PathBuf::from("data/ref").join(fasta_dist);

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
    fasta_path
}

/// Download a bam file if it's an url, otherwise check the file exist
pub fn resolve_bam(bam: &str, outdir: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    if bam.starts_with("http://") || bam.starts_with("https://") {
        download_bam(bam, outdir)
    } else {
        let path = PathBuf::from(bam);
        if !path.exists() {
            return Err(format!("BAM file not found: {}", path.display()).into());
        }
        Ok(path)
    }
}

fn download_bam(url: &str, outdir: &PathBuf) -> Result<PathBuf, Box<dyn Error>> {
    let parsed = Url::parse(&url)?;
    log::info!("Downloading bam file {}", url);
    let filename = parsed
        .path_segments()
        .and_then(|s| s.last())
        .filter(|s| !s.is_empty())
        .ok_or("could not extract filename from URL")?;

    let output = outdir.join(filename);
    download_blocking(&url, &output);
    Ok(output)
}
/// Generate controls from clinvar variant into a BAM file
/// 1. Sample clinvar randomly
/// 2. Generate mutation files from those clinvar variant for varben
/// 3. Apply varben on a single bam file
///
/// If `clinvar_vcf` is None the file is downloaded from NCBI and stored in `outdir`.
pub fn edit_bam(
    bam: &PathBuf,
    bed: &PathBuf,
    capture: &str,
    fasta: PathBuf,
    clinvar_vcf: Option<PathBuf>,
    nb_variants: Option<u32>,
    outdir: PathBuf,
) -> Result<PathBuf, Box<dyn Error>> {
    std::fs::create_dir_all(&outdir)?;

    index_bam(bam)?;
    let bam_cleaned = remove_hard_clips(bam)?;

    let clinvar_path = match clinvar_vcf {
        Some(p) => p,
        None => download_clinvar(),
    };

    let n = match nb_variants {
        Some(n) => n,
        None => 1000,
    };

    // sample n clinvar variants 50bp apart
    let vcf_out = outdir.join(format!("clinvar_{capture}.vcf.gz"));
    let mut_out = outdir.join(format!("clinvar_{capture}.mut"));
    let header = sample_clinvar(
        clinvar_path,
        PathBuf::from(bed),
        50,
        n,
        vcf_out,
        mut_out.clone(),
    )?;
    let edited_bam = insert_variants(mut_out, bam_cleaned, outdir.clone(), &fasta)?;

    write_varben_as_vcf(bam, header, outdir, &fasta)?;

    Ok(edited_bam)
}

/// Use varben (muteditor) to insert a list of variant in a bed files. Require varben, samtools, bwa
/// Require a reference genome and a bwa index
/// Output folder is the varben subfolder of `outdir`
fn insert_variants(
    mut_file: PathBuf,
    bam: PathBuf,
    outdir: PathBuf,
    fasta: &PathBuf,
) -> Result<PathBuf, Box<dyn Error>> {
    let outdir_ = outdir.join("varben");
    std::fs::create_dir_all(&outdir_)?;

    let fasta_str = fasta.to_str().ok_or("Invalid fasta path")?;
    let mut_str = mut_file.to_str().ok_or("Invalid mutation file path")?;
    let bam_str = bam.to_str().ok_or("Invalid BAM path")?;
    let outdir_str = outdir_.to_str().ok_or("Invalid output folder")?;

    let output = outdir_.join("edit.sorted.bam");
    if output.exists() {
        log::debug!("{:?} already exists", output);
    } else {
        let log_path = outdir_.join("varben.log");
        let status = insert_variants_varben(mut_str, bam_str, fasta_str, outdir_str, log_path)?;
        if !status.success() {
            return Err(format!("muteditor exited with status {status}").into());
        }
    }
    Ok(output)
}

fn insert_variants_varben(
    mut_str: &str,
    bam_str: &str,
    fasta_str: &str,
    outdir_str: &str,
    log_path: PathBuf,
) -> Result<ExitStatus, std::io::Error> {
    let log_file = File::create(&log_path)?;

    Command::new("muteditor")
        .args([
            "-m",
            mut_str,
            "-p",
            &nb_threads().to_string(),
            "-b",
            bam_str,
            "-r",
            fasta_str,
            "--aligner",
            "bwa",
            "--alignerIndex",
            fasta_str,
            "-o",
            outdir_str,
        ])
        .stdout(Stdio::from(log_file))
        .status()
}
/// samtools fastq is the best tool in our testing. It requires read to be sorted by read name beforehand
/// Fastq output will be in the same directory as the BAM and suffixed with _1.fq.gz
pub fn bam_to_fastq(
    bam: PathBuf,
    fq1: PathBuf,
    fq2: PathBuf,
) -> Result<(PathBuf, PathBuf), Box<dyn Error>> {
    let cmd = format!(
        "samtools sort -n {bam} -@ {threads} | samtools fastq -1 {fq1} -2 {fq2} -0 /dev/null -s /dev/null -n -@ {threads} -",
        bam = bam.display(),
        fq1 = fq1.display(),
        fq2 = fq2.display(),
        threads = nb_threads(),
    );

    let status = Command::new("sh").args(["-c", &cmd]).status()?;

    if !status.success() {
        return Err(format!("bam to fastq conversion exited with {status}").into());
    }
    Ok((fq1, fq2))
}

fn nb_threads() -> usize {
    thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

/// Convert varben tabular data to a simple VCF.
/// Require a fasta as varben only store alternative allel and not refercen
fn write_varben_as_vcf_single(
    mut_file: PathBuf,
    vcf_out: PathBuf,
    header: vcf::Header,
    fasta: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let mut fasta_reader = fasta::io::indexed_reader::Builder::default().build_from_path(fasta)?;

    let reader = BufReader::new(File::open(mut_file)?);
    let mut records: Vec<RecordBuf> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 6 {
            continue;
        }
        let chrom = fields[0];
        let pos: usize = fields[1].parse()?;
        let alt = fields[5].trim();

        let region: noodles::core::Region = format!("{}:{}-{}", chrom, pos, pos).parse()?;
        let record = fasta_reader.query(&region)?;
        let seq = record.sequence().as_ref();
        let ref_base = std::str::from_utf8(&seq[..1])?.to_uppercase();

        let buf = RecordBuf::builder()
            .set_reference_sequence_name(chrom)
            .set_variant_start(noodles::core::Position::try_from(pos)?)
            .set_reference_bases(ref_base)
            .set_alternate_bases(AltBasesBuf::from(vec![alt.to_string()]))
            .build();

        records.push(buf);
    }

    sort_by_chromosome(&mut records);

    // let header = vcf::Header::default();
    let mut writer = File::create(vcf_out.clone())
        .map(bgzf::io::Writer::new)
        .map(vcf::io::Writer::new)?;
    writer.write_header(&header)?;
    for rec in &records {
        writer.write_variant_record(&header, rec)?;
    }

    log::debug!("Wrote {} variants to {:?}", records.len(), vcf_out);
    Ok(())
}
