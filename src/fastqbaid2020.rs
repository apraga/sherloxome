//! # Raw FASTQ data according to baid 2020
//!
//! Different combinations are available acoording to GIAB patients,  sequencer types, capture kit
//! Note : not all combination are availabe in Baid's data (for depth mostly).
//! For GIAB patients, [see here](crate::giab)
//!
use serde::Deserialize;

/// Capture kit
/// - Agilent_SureSelect_All_Exons_v7_hg38
/// - Truseq exome: we use the version lifter in hg38 through UCSC
/// - IDT-xGen : xgen-exome-hyb-panel-v2-targets-hg38
#[derive(Deserialize, Debug)]
pub enum Kit {
    Agilent,
    Idt,
    Truseq,
}

/// Sequencer (novaseq, hiseq)
#[derive(Deserialize, Debug)]
pub enum Sequencer {
    Novaseq,
    Hiseq400,
}

/// Depth cannot be defined by a number, so the `DP` prefix is added
#[derive(Deserialize, Debug)]
pub enum Depth {
    DP50,
    DP75,
    DP100,
}
