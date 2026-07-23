use crate::silico::nb_threads;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::variant::record::AlternateBases;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

/// Write one SNV in simuscop tab-separated variation format (no header line).
/// Always emits `het` since we sample heterozygous AF from clinvar.
fn write_variant(
    w: &mut impl Write,
    record: &RecordBuf,
    sample: &str,
) -> Result<(), Box<dyn Error>> {
    let chrom = record.reference_sequence_name();
    let pos = usize::from(record.variant_start().unwrap());
    let ref_base = record.reference_bases().to_lowercase();
    let alt = record
        .alternate_bases()
        .iter()
        .filter_map(|a| a.ok())
        .next()
        .unwrap_or(".");
    writeln!(w, "s\t{sample}\t{chrom}\t{pos}\t{ref_base}\t{alt}\thet")?;
    Ok(())
}

/// Write simuscop input file from a list of variant
pub fn write_input(variants: &[RecordBuf], out: &Path, sample: &str) -> Result<(), Box<dyn Error>> {
    let mut w = BufWriter::new(File::create(out)?);
    for v in variants {
        write_variant(&mut w, v, sample)?;
    }
    w.flush()?;
    Ok(())
}

/// Write a simuReads config file
pub fn write_config(
    config_path: &Path,
    fasta: &Path,
    profile_dir: &Path,
    variation: &Path,
    bed: &Path,
    name: &str,
    output_dir: &Path,
    coverage: u32,
) -> Result<(), Box<dyn Error>> {
    let threads = nb_threads();
    let mut w = BufWriter::new(File::create(config_path)?);
    writeln!(w, "ref = {}", fasta.display())?;
    writeln!(w, "profile = {}", profile_dir.display())?;
    writeln!(w, "variation = {}", variation.display())?;
    writeln!(w, "target = {}", bed.display())?;
    writeln!(w, "name = {}", name)?;
    writeln!(w, "output = {}", output_dir.display())?;
    writeln!(w, "layout = PE")?;
    writeln!(w, "threads = {}", threads)?;
    writeln!(w, "coverage = {}", coverage)?;
    w.flush()?;
    Ok(())
}
