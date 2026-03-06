//! # GIAB patients
//!
//! Here are defined the patients available in GIAB reference dataset. This module provide access to the VCF.
//! For raw data, [see here](crate::fastqbaid2020)
use serde::Deserialize;
use std::fmt;

/// Patient ID according to GIAB
/// Equality is required to compare runs [see here](crate::fastqbaid2020)
#[derive(Clone, Copy, Deserialize, Debug, Hash, Eq, PartialEq)]
pub enum Patient {
    HG001,
    HG002,
    HG003,
    HG004,
    HG005,
    HG006,
    HG007,
}

impl fmt::Display for Patient {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Patient::HG001 => "HG001",
            Patient::HG002 => "HG002",
            Patient::HG003 => "HG003",
            Patient::HG004 => "HG004",
            Patient::HG005 => "HG005",
            Patient::HG006 => "HG006",
            Patient::HG007 => "HG007",
        };
        write!(f, "{s}")
    }
}

const BASE_URL: &str = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release";

/// GIAB FTP sub-path that is unique to each patient.
fn url_tail(p: &Patient) -> &'static str {
    match p {
        Patient::HG001 => "NA12878_HG001/latest/GRCh38",
        Patient::HG002 => "AshkenazimTrio/HG002_NA24385_son/latest/GRCh38",
        Patient::HG003 => "AshkenazimTrio/HG003_NA24149_father/latest/GRCh38",
        Patient::HG004 => "AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38",
        Patient::HG005 => "ChineseTrio/HG005_NA24631_son/latest/GRCh38",
        Patient::HG006 => "ChineseTrio/HG006_NA24694_father/latest/GRCh38",
        Patient::HG007 => "ChineseTrio/HG007_NA24695_mother/latest/GRCh38",
    }
}

/// GIAB full URL on their FTP
pub fn vcf_url(p: &Patient) -> String {
    format!("{BASE_URL}/{}/{}", url_tail(p), vcf_file(p))
}

/// GIAB file name on their FTP
pub fn vcf_file(p: &Patient) -> String {
    format!("{p}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
}

/// GIAB Base file name for the benchmark BED.
/// HG001–HG004 ship a `_noinconsistent` BED; HG005–HG007 use the plain one.
pub fn bed_file(p: &Patient) -> String {
    let suffix = match p {
        Patient::HG002 | Patient::HG003 | Patient::HG004 => "_noinconsistent.bed",
        _ => ".bed",
    };
    format!("{p}_GRCh38_1_22_v4.2.1_benchmark{suffix}")
}

/// GIAB full URL on their FTP
pub fn bed_url(p: &Patient) -> String {
    format!("{BASE_URL}/{}/{}", url_tail(p), bed_file(p))
}
