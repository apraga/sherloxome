//! # Raw FASTQ data according to baid 2020
//!
//! Different combinations are available acoording to GIAB patients,  sequencer types, capture kit
//! Note : not all combination are availabe in Baid's data (for depth mostly).
//! For GIAB patients, [see here](crate::giab)
//!
use crate::giab::{Patient, patient_from_filename};
use crate::setup::SamplesheetRow;
use serde::Deserialize;
use std::path::PathBuf;
use std::{collections::HashSet, fmt};

/// A run is defined by a [Patient], [Sequencer], capture [Capture] and [Depth]
#[derive(Clone, Copy, Deserialize, Debug, Hash, Eq, PartialEq)]
#[serde(rename_all = "lowercase")]
pub struct Run {
    pub patient: Patient,
    pub capture: Capture,
    pub sequencer: Sequencer,
    pub depth: Depth,
}

/// Capture kit
/// - Agilent_SureSelect_All_Exons_v7_hg38
/// - Truseq exome: we use the version lifter in hg38 through UCSC
/// - IDT-xGen : xgen-exome-hyb-panel-v2-targets-hg38
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

pub fn all_captures() -> Vec<Capture> {
    [Capture::Agilent, Capture::Idt, Capture::Truseq].to_vec()
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

pub fn all_sequencers() -> Vec<Sequencer> {
    [Sequencer::Hiseq4000, Sequencer::Novaseq].to_vec()
}

/// Depth cannot be defined by a number, so the `DP` prefix is added
#[derive(Copy, Clone, Deserialize, Debug, Hash, Eq, PartialEq)]
pub enum Depth {
    #[serde(rename = "50x")]
    DP50,
    #[serde(rename = "75x")]
    DP75,
    #[serde(rename = "100x")]
    DP100,
}

impl fmt::Display for Depth {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Depth::DP50 => "50x",
            Depth::DP75 => "75x",
            Depth::DP100 => "100x",
        };
        write!(f, "{s}")
    }
}
pub fn all_depths() -> Vec<Depth> {
    [Depth::DP50, Depth::DP75, Depth::DP100].to_vec()
}

/// Return all available combinations for filtering later on.
/// We have to hardcode it as there is no simple rules
pub fn available() -> HashSet<Run> {
    let combinations = [
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Hiseq4000,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            capture: Capture::Truseq,
            depth: Depth::DP75,
        },
    ];

    HashSet::from(combinations)
}

/// All capture captureavailable in Baid2020
pub fn capture() -> Vec<Capture> {
    [Capture::Agilent, Capture::Idt, Capture::Truseq].to_vec()
}

pub fn real_row(run: &Run) -> SamplesheetRow {
    let sample = format!(
        "{:}-{:}-{:}-{:}",
        run.patient, run.sequencer, run.capture, run.depth
    );

    return SamplesheetRow {
        patient: run.patient.to_string(),
        sample: sample,
        lane: 1,
        fastq_1: url(*run, "R1"),
        fastq_2: url(*run, "R2"),
    };
}

/// Use google cloud URL. Nextflow will download the data
pub fn url(run: Run, lane: &str) -> String {
    let root = format!(
        "https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/{sequencer}/wes_{capture}/{depth}/{patient}.{sequencer}.wes_{capture}.{depth}",
        capture = run.capture,
        sequencer = run.sequencer,
        patient = run.patient,
        depth = run.depth,
    );

    format!("{}.{}.fastq.gz", root, lane)
}

pub fn run_from_filename(fname: &PathBuf) -> Option<Run> {
    let name = fname.to_string_lossy();
    let patient = patient_from_filename(fname);
    let capture = all_captures()
        .into_iter()
        .find(|p| name.contains(&p.to_string()));
    let depth = all_depths()
        .into_iter()
        .find(|p| name.contains(&p.to_string()));
    let sequencer = all_sequencers()
        .into_iter()
        .find(|p| name.contains(&p.to_string()));
    Some(Run {
        patient: patient?,
        capture: capture?,
        sequencer: sequencer?,
        depth: depth?,
    })
}
