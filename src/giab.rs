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
