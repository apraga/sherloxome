//! # Insilico Controls
//! Generate control variants by sampling clinvar data

use crate::download_blocking;
use log;
use noodles::bed;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::info::field::Value;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

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

/// Select all clinvar pathogenic SNV inside the capture kit 50bp apart
/// Outputpath is hardcoded to data/exp_raw/clinvar_$CAPTURE
/// The more efficient way is to parse the VCF and check if each variant is in the BED (stored as an interval tree)
pub fn sample_clinvar(bed: PathBuf, spacing: u32, output: PathBuf) -> Result<(), Box<dyn Error>> {
    let vcf = download_clinvar();

    // Load BED into interval tree per chrom
    let capture = load_bed(&bed)?;
    let mut vcf_reader = File::open(vcf)
        .map(bgzf::io::Reader::new)
        .map(vcf::io::Reader::new)?;
    let vcf_header = vcf_reader.read_header()?;
    let mut vcf_writer = File::create(&output)
        .map(bgzf::io::Writer::new)
        .map(vcf::io::Writer::new)?;
    vcf_writer.write_header(&vcf_header)?;

    let mut last_pos: HashMap<String, u64> = HashMap::new();

    for result in vcf_reader.records() {
        let record = result?;
        if let Some((chrom, pos)) = keep_variant(&record, &vcf_header, &capture, &last_pos, spacing)
        {
            last_pos.insert(chrom, pos);
            vcf_writer.write_record(&vcf_header, &record)?;
        }
    }
    log::debug!("Wrote {:?}", output);
    Ok(())
}

/// Do it for a single record

fn keep_variant(
    record: &vcf::Record,
    header: &vcf::Header,
    capture: &BedIndex,
    last_pos: &HashMap<String, u64>,
    spacing: u32,
) -> Option<(String, u64)> {
    let chrom = record.reference_sequence_name().to_string();
    let pos = usize::from(record.variant_start()?.ok()?) as u64;
    let last = last_pos.get(&chrom).copied().unwrap_or(0);

    (in_capture(capture, &chrom, pos)
        && is_not_benign(&record.info(), header)
        && is_snv(record)
        && pos != last
        && pos.abs_diff(last) >= spacing.into())
    .then_some((chrom, pos))
}
