//! # GIAB patients
//!
//! Here are defined the patients available in GIAB reference dataset. This module provide access to the VCF.
//! For raw data, [see here](crate::fastqbaid2020)
use serde::Deserialize;

/// Patient ID
#[derive(Deserialize, Debug)]
pub enum Patient {
    HG001,
    HG002,
    HG003,
    HG004,
    HG005,
    HG006,
    HG007,
}
