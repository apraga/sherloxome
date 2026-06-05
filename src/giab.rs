//! # GIAB patients
//!
//! Here are defined the patients available in GIAB reference dataset. This module provide access to the VCF.
//! For raw data, [see here](crate::fastqbaid2020)
use crate::ref_dir;
use serde::Deserialize;
use std::fmt;
use std::path::PathBuf;
use std::str::FromStr;

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

impl FromStr for Patient {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "HG001" => Ok(Patient::HG001),
            "HG002" => Ok(Patient::HG002),
            "HG003" => Ok(Patient::HG003),
            "HG004" => Ok(Patient::HG004),
            "HG005" => Ok(Patient::HG005),
            "HG006" => Ok(Patient::HG006),
            "HG007" => Ok(Patient::HG007),
            _ => Err(()),
        }
    }
}

const BASE_URL: &str = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release";

/// GIAB FTP sub-path that is unique to each patient.
fn url_tail(p: &Patient) -> &'static str {
    match p {
        Patient::HG001 => "NA12878_HG001/NISTv4.2.1/GRCh38",
        Patient::HG002 => "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38",
        Patient::HG003 => "AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38",
        Patient::HG004 => "AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38",
        Patient::HG005 => "ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38",
        Patient::HG006 => "ChineseTrio/HG006_NA24694_father/NISTv4.2.1/GRCh38",
        Patient::HG007 => "ChineseTrio/HG007_NA24695_mother/NISTv4.2.1/GRCh38",
    }
}

/// GIAB full URL on their FTP
pub fn vcf_url(p: &Patient) -> String {
    format!("{BASE_URL}/{}/{}", url_tail(p), vcf_filename(p).display())
}

/// GIAB file name on their FTP
pub fn vcf_filename(p: &Patient) -> PathBuf {
    PathBuf::from(format!("{p}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"))
}

/// Full path locally (in data/ref). Must be download beforehand
pub fn vcf_path(p: &Patient) -> PathBuf {
    ref_dir().join(vcf_filename(p))
}

/// GIAB Base file name for the benchmark BED.
/// HG001–HG004 ship a `_noinconsistent` BED; HG005–HG007 use the plain one.
pub fn bed_filename(p: &Patient) -> PathBuf {
    let suffix = match p {
        Patient::HG002 | Patient::HG003 | Patient::HG004 => "_noinconsistent.bed",
        _ => ".bed",
    };
    PathBuf::from(format!("{p}_GRCh38_1_22_v4.2.1_benchmark{suffix}"))
}

/// Full GIAB bed file locally
pub fn bed_path(p: &Patient) -> PathBuf {
    ref_dir().join(bed_filename(p))
}

/// GIAB full URL on their FTP
pub fn bed_url(p: &Patient) -> String {
    format!("{BASE_URL}/{}/{}", url_tail(p), bed_filename(p).display())
}

/// GIAB full URL for the VCF tabix index (.tbi)
pub fn tbi_url(p: &Patient) -> String {
    format!("{}.tbi", vcf_url(p))
}

/// Full path locally for the VCF tabix index
pub fn tbi_path(p: &Patient) -> PathBuf {
    PathBuf::from(format!("{}.tbi", vcf_path(p).display()))
}

/// Infer the GIAB [`Patient`] from a file path by matching its string representation.
pub fn patient_from_filename(fname: &PathBuf) -> Option<Patient> {
    let name = fname.to_string_lossy();
    all_patients()
        .into_iter()
        .find(|p| name.contains(&p.to_string()))
}

/// All GIAB patients
pub fn all_patients() -> Vec<Patient> {
    [
        Patient::HG001,
        Patient::HG002,
        Patient::HG003,
        Patient::HG004,
        Patient::HG005,
        Patient::HG006,
        Patient::HG007,
    ]
    .to_vec()
}
