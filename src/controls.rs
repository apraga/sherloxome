//! # Insilico Controls
//! Generate control variants by sampling clinvar data

use crate::download_blocking;
use log;
use noodles::bed;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::info::field::Value;
use rand::RngExt;
use rand::prelude::IteratorRandom;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;

struct Snv {
    chrom: String, // chr-prefixed
    pos: u64,      // 1-based
    alt: String,
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
        let chrom = std::str::from_utf8(record.reference_sequence_name())?
            .trim_start_matches("chr")
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

// Generate AF betwen 0.4 and 0.6
fn write_varben(
    w: &mut impl Write,
    snv: &Snv,
    rng: &mut impl RngExt,
) -> Result<(), Box<dyn Error>> {
    let af: f64 = rng.random_range(0.4..0.6);
    writeln!(
        w,
        "{}\t{}\t{}\t{}\tsnv\t{}",
        snv.chrom, snv.pos, snv.pos, af, snv.alt
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

fn sample_clinvar_variants(
    reader: &mut vcf::io::Reader<bgzf::io::Reader<File>>,
    header: &vcf::Header,
    capture: &BedIndex,
    spacing: u32,
    n: usize,
    rng: &mut impl RngExt,
) -> Vec<(vcf::Record, Snv)> {
    let mut last_pos: HashMap<String, u64> = HashMap::new();
    let mut selected: Vec<(vcf::Record, Snv)> = reader
        .records()
        .filter_map(|r| r.ok())
        .filter_map(|record| {
            keep_variant(&record, &header, &capture, &last_pos, spacing).map(|snv| {
                last_pos.insert(snv.chrom.clone(), snv.pos);
                (record, snv)
            })
        })
        .sample(rng, n);
    selected.sort_by(|(_, a), (_, b)| {
        chrom_order(&a.chrom)
            .cmp(&chrom_order(&b.chrom))
            .then(a.pos.cmp(&b.pos))
    });

    log::debug!("Selected {} variants", selected.len());
    selected
}

/// Select n clinvar pathogenic SNV inside the capture kit 50bp apart.
/// Output is both a mutation file for varben and a VCF (hardocoded into data/exp_raw/clinvar_$CAPTURE)
/// The most efficient way is to parse clinvar once to get variant in the capture kit and 50pb apart.
/// In a second pass, sample randomly n of theme
pub fn sample_clinvar(
    bed: PathBuf,
    spacing: u32,
    n: usize,
    vcf_out: PathBuf,
    mut_out: PathBuf,
) -> Result<(), Box<dyn Error>> {
    let vcf = download_clinvar();

    // Prepare VCF file
    // Load BED into interval tree per chrom
    let capture = load_bed(&bed)?;
    let mut reader = File::open(vcf)
        .map(bgzf::io::Reader::new)
        .map(vcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut rng = rand::rng();
    let variants = sample_clinvar_variants(&mut reader, &header, &capture, spacing, n, &mut rng);

    write_sampled_clinvar_vcf(&variants, &header, &vcf_out)?;
    write_sampled_clinvar_mut(&variants, &mut_out, &mut rng)?;
    log::debug!("Wrote {:?} and {:?}", vcf_out, mut_out);
    Ok(())
}

fn write_sampled_clinvar_vcf(
    variants: &Vec<(vcf::Record, Snv)>,
    header: &vcf::Header,
    vcf_out: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let mut writer = File::create(&vcf_out)
        .map(bgzf::io::Writer::new)
        .map(vcf::io::Writer::new)?;
    writer.write_header(&header)?;

    for (record, _) in variants {
        writer.write_record(&header, record)?;
    }
    Ok(())
}

fn write_sampled_clinvar_mut(
    variants: &Vec<(vcf::Record, Snv)>,
    mut_out: &PathBuf,
    rng: &mut impl RngExt,
) -> Result<(), Box<dyn Error>> {
    // Prepare mutation file
    let mut varben_w = BufWriter::new(File::create(&mut_out)?);
    writeln!(varben_w, "#chrom\tstart\tend\tAF\ttype\talt")?;

    for (_, snv) in variants {
        write_varben(&mut varben_w, snv, rng)?;
    }
    varben_w.flush()?;
    Ok(())
}
/// Do it for a single record
/// Input VCF may use chr prefix, output will have `chr` prefix.
fn keep_variant(
    record: &vcf::Record,
    header: &vcf::Header,
    capture: &BedIndex,
    last_pos: &HashMap<String, u64>,
    spacing: u32,
) -> Option<Snv> {
    let chrom_raw = record.reference_sequence_name().to_string();
    // Remove chr prefix if needed
    let chrom_nb = chrom_raw.strip_prefix("chr").unwrap_or(&chrom_raw);
    let chrom_chr = format!("chr{}", chrom_nb);

    let pos = usize::from(record.variant_start()?.ok()?) as u64;
    let last = last_pos.get(chrom_nb).copied().unwrap_or(0);
    let alt = record
        .alternate_bases()
        .iter()
        .filter_map(|a| a.ok())
        .next()
        .unwrap_or(".")
        .to_string();
    (in_capture(capture, &chrom_nb, pos)
        && is_not_benign(&record.info(), header)
        && is_snv(record)
        && pos != last
        && pos.abs_diff(last) >= spacing.into())
    .then_some(Snv {
        chrom: chrom_chr,
        pos,
        alt,
    })
}
