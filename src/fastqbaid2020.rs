//! # Raw FASTQ data according to baid 2020
//!
//! Different combinations are available acoording to GIAB patients,  sequencer types, capture kit
//! Note : not all combination are availabe in Baid's data (for depth mostly).
//! For GIAB patients, [see here](crate::giab)
//!
use crate::giab::Patient;
use serde::Deserialize;
use std::fs::File;
use std::io::Write;
use std::{collections::HashSet, fmt};

/// A run is defined by a [Patient], [Sequencer], capture [Kit] and [Depth]
#[derive(Clone, Copy, Deserialize, Debug, Hash, Eq, PartialEq)]
#[serde(rename_all = "lowercase")]
pub struct Run {
    pub patient: Patient,
    pub kit: Kit,
    pub sequencer: Sequencer,
    pub depth: Depth,
}

/// Capture kit
/// - Agilent_SureSelect_All_Exons_v7_hg38
/// - Truseq exome: we use the version lifter in hg38 through UCSC
/// - IDT-xGen : xgen-exome-hyb-panel-v2-targets-hg38
#[derive(Copy, Clone, Deserialize, Debug, Hash, Eq, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Kit {
    Agilent,
    Idt,
    Truseq,
}

impl fmt::Display for Kit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Kit::Agilent => "agilent",
            Kit::Idt => "idt",
            Kit::Truseq => "truseq",
        };
        write!(f, "{s}")
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

/// Return all available combinations for filtering later on.
/// We have to hardcode it as there is no simple rules
pub fn available() -> HashSet<Run> {
    let combinations = [
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG001,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG002,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG003,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG004,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG005,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG006,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Hiseq4000,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Agilent,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP100,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Idt,
            depth: Depth::DP75,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP50,
        },
        Run {
            patient: Patient::HG007,
            sequencer: Sequencer::Novaseq,
            kit: Kit::Truseq,
            depth: Depth::DP75,
        },
    ];

    HashSet::from(combinations)
}

pub fn samplesheet_real(runs: HashSet<Run>) {
    let mut file = File::create("samplesheet.csv").unwrap();
    file.write(b"patient,sample,lane,fastq_1,fastq_2\n")
        .unwrap();
    for run in runs {
        samplesheet_row(run, &mut file);
    }
}

pub fn samplesheet_row(run: Run, file: &mut File) {
    let sample = format!(
        "{:}-{:}-{:}-{:}",
        run.patient, run.sequencer, run.kit, run.depth
    );
    println!("{:?}", sample);
    let fastq1 = url(run, "R1");
    let fastq2 = url(run, "R2");

    let row = format!("{},{sample},1,{fastq1},{fastq2}\n", run.patient);
    file.write(row.as_bytes()).unwrap();
}

/// Use google cloud URL. Nextflow will download the data
pub fn url(run: Run, lane: &str) -> String {
    let root = format!(
        "https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/{sequencer}/wes_{kit}/{depth}/{patient}.{sequencer}.wes_{kit}.{depth}",
        kit = run.kit,
        sequencer = run.sequencer,
        patient = run.patient,
        depth = run.depth,
    );

    format!("{}.{}.fastq.gz", root, lane)
}
