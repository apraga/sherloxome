//! # Raw FASTQ data according to baid 2020
//!
//! Different combinations are available acoording to GIAB patients,  sequencer types, capture kit. `crate::run` defines a filenaming scheme but without enforcing a set of values.
//! Here the number of values is limited to the paper data.
//! Note : not all combination are availabe in Baid's data (for depth mostly).
//! For GIAB patients, [see here](crate::giab)
//!
use crate::giab::Patient;
use crate::run::Run;
use crate::run::run_to_string;
use crate::setup::SamplesheetRow;
use serde::Deserialize;
use std::collections::HashSet;
use std::fmt;
use std::str::FromStr;

/// Capture kit
/// - Agilent_SureSelect_All_Exons_v7_hg38
/// - Agilent limited to COL6A1 (for testing only)
/// - Truseq exome: we use the version lifter in hg38 through UCSC
/// - IDT-xGen : xgen-exome-hyb-panel-v2-targets-hg38
#[derive(Copy, Clone, Deserialize, Debug, Hash, Eq, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Capture {
    Agilent,
    Idt,
    Truseq,
}

impl fmt::Display for Capture {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Capture::Agilent => "agilent",
            Capture::Idt => "idt",
            Capture::Truseq => "truseq",
        };
        write!(f, "{s}")
    }
}

impl FromStr for Capture {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "agilent" => Ok(Capture::Agilent),
            "idt" => Ok(Capture::Idt),
            "truseq" => Ok(Capture::Truseq),
            _ => Err(()),
        }
    }
}

/// Sequencer (novaseq, hiseq)
#[derive(Copy, Clone, Deserialize, Debug, Hash, Eq, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Sequencer {
    Hiseq4000,
    Novaseq,
}

impl fmt::Display for Sequencer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Sequencer::Hiseq4000 => "hiseq4000",
            Sequencer::Novaseq => "novaseq",
        };
        write!(f, "{s}")
    }
}

impl FromStr for Sequencer {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "hiseq4000" => Ok(Sequencer::Hiseq4000),
            "novaseq" => Ok(Sequencer::Novaseq),
            _ => Err(()),
        }
    }
}

/// All supported sequencing depths.
pub fn all_depths() -> Vec<u32> {
    [50, 75, 100].to_vec()
}

/// Convert Baid202 types to string representation of `crate::run`
pub fn to_run(p: Patient, s: Sequencer, c: Capture, d: u32, silico: Option<String>) -> Run {
    Run {
        sample: format!("{}", p),
        sequencer: format!("{}", s),
        capture: format!("{}", c),
        depth: d,
        silico: silico,
    }
}
/// Return all available combinations for filtering later on.
/// We have to hardcode it as there is no simple rules
pub fn available() -> HashSet<Run> {
    let combinations = [
        to_run(
            Patient::HG001,
            Sequencer::Hiseq4000,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(Patient::HG001, Sequencer::Hiseq4000, Capture::Idt, 50, None),
        to_run(
            Patient::HG001,
            Sequencer::Hiseq4000,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG001,
            Sequencer::Novaseq,
            Capture::Agilent,
            100,
            None,
        ),
        to_run(
            Patient::HG001,
            Sequencer::Novaseq,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(
            Patient::HG001,
            Sequencer::Novaseq,
            Capture::Agilent,
            75,
            None,
        ),
        to_run(
            Patient::HG001,
            Sequencer::Novaseq,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG001,
            Sequencer::Novaseq,
            Capture::Truseq,
            75,
            None,
        ),
        to_run(
            Patient::HG002,
            Sequencer::Hiseq4000,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(Patient::HG002, Sequencer::Hiseq4000, Capture::Idt, 50, None),
        to_run(
            Patient::HG002,
            Sequencer::Hiseq4000,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG002,
            Sequencer::Novaseq,
            Capture::Agilent,
            100,
            None,
        ),
        to_run(
            Patient::HG002,
            Sequencer::Novaseq,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(
            Patient::HG002,
            Sequencer::Novaseq,
            Capture::Agilent,
            75,
            None,
        ),
        to_run(
            Patient::HG002,
            Sequencer::Novaseq,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG002,
            Sequencer::Novaseq,
            Capture::Truseq,
            75,
            None,
        ),
        to_run(
            Patient::HG003,
            Sequencer::Hiseq4000,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(Patient::HG003, Sequencer::Hiseq4000, Capture::Idt, 50, None),
        to_run(
            Patient::HG003,
            Sequencer::Hiseq4000,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG003,
            Sequencer::Novaseq,
            Capture::Agilent,
            100,
            None,
        ),
        to_run(
            Patient::HG003,
            Sequencer::Novaseq,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(
            Patient::HG003,
            Sequencer::Novaseq,
            Capture::Agilent,
            75,
            None,
        ),
        to_run(
            Patient::HG003,
            Sequencer::Novaseq,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG003,
            Sequencer::Novaseq,
            Capture::Truseq,
            75,
            None,
        ),
        to_run(
            Patient::HG004,
            Sequencer::Hiseq4000,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(Patient::HG004, Sequencer::Hiseq4000, Capture::Idt, 50, None),
        to_run(
            Patient::HG004,
            Sequencer::Hiseq4000,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG004,
            Sequencer::Novaseq,
            Capture::Agilent,
            100,
            None,
        ),
        to_run(
            Patient::HG004,
            Sequencer::Novaseq,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(
            Patient::HG004,
            Sequencer::Novaseq,
            Capture::Agilent,
            75,
            None,
        ),
        to_run(
            Patient::HG004,
            Sequencer::Novaseq,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG004,
            Sequencer::Novaseq,
            Capture::Truseq,
            75,
            None,
        ),
        to_run(
            Patient::HG005,
            Sequencer::Hiseq4000,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(Patient::HG005, Sequencer::Hiseq4000, Capture::Idt, 50, None),
        to_run(
            Patient::HG005,
            Sequencer::Hiseq4000,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG005,
            Sequencer::Novaseq,
            Capture::Agilent,
            100,
            None,
        ),
        to_run(
            Patient::HG005,
            Sequencer::Novaseq,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(
            Patient::HG005,
            Sequencer::Novaseq,
            Capture::Agilent,
            75,
            None,
        ),
        to_run(
            Patient::HG005,
            Sequencer::Novaseq,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG005,
            Sequencer::Novaseq,
            Capture::Truseq,
            75,
            None,
        ),
        to_run(
            Patient::HG006,
            Sequencer::Hiseq4000,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(Patient::HG006, Sequencer::Hiseq4000, Capture::Idt, 50, None),
        to_run(
            Patient::HG006,
            Sequencer::Hiseq4000,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG006,
            Sequencer::Novaseq,
            Capture::Agilent,
            100,
            None,
        ),
        to_run(
            Patient::HG006,
            Sequencer::Novaseq,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(
            Patient::HG006,
            Sequencer::Novaseq,
            Capture::Agilent,
            75,
            None,
        ),
        to_run(
            Patient::HG006,
            Sequencer::Novaseq,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG006,
            Sequencer::Novaseq,
            Capture::Truseq,
            75,
            None,
        ),
        to_run(
            Patient::HG007,
            Sequencer::Hiseq4000,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(Patient::HG007, Sequencer::Hiseq4000, Capture::Idt, 50, None),
        to_run(
            Patient::HG007,
            Sequencer::Hiseq4000,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG007,
            Sequencer::Novaseq,
            Capture::Agilent,
            100,
            None,
        ),
        to_run(
            Patient::HG007,
            Sequencer::Novaseq,
            Capture::Agilent,
            50,
            None,
        ),
        to_run(
            Patient::HG007,
            Sequencer::Novaseq,
            Capture::Agilent,
            75,
            None,
        ),
        to_run(
            Patient::HG007,
            Sequencer::Novaseq,
            Capture::Truseq,
            50,
            None,
        ),
        to_run(
            Patient::HG007,
            Sequencer::Novaseq,
            Capture::Truseq,
            75,
            None,
        ),
    ];

    HashSet::from(combinations)
}

/// All capture captureavailable in Baid2020
pub fn capture() -> Vec<Capture> {
    [Capture::Agilent, Capture::Idt, Capture::Truseq].to_vec()
}

pub fn real_row(run: &Run) -> SamplesheetRow {
    return SamplesheetRow {
        patient: run.sample.to_string(),
        sample: run_to_string(run),
        lane: 1,
        fastq_1: url(run, "R1"),
        fastq_2: url(run, "R2"),
    };
}

/// Use google cloud URL. Nextflow will download the data
pub fn url(run: &Run, lane: &str) -> String {
    let root = format!(
        "https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/{sequencer}/wes_{capture}/{depth}/{sample}.{sequencer}.wes_{capture}.{depth}",
        capture = run.capture,
        sequencer = run.sequencer,
        sample = run.sample,
        depth = run.depth,
    );

    format!("{}.{}.fastq.gz", root, lane)
}

/// All supported capture kits.
pub fn all_captures() -> Vec<Capture> {
    [
        Capture::Agilent,
        Capture::AgilentCOL6A1,
        Capture::Idt,
        Capture::Truseq,
    ]
    .to_vec()
}

/// All supported sequencer types.
pub fn all_sequencers() -> Vec<Sequencer> {
    [Sequencer::Hiseq4000, Sequencer::Novaseq].to_vec()
}
