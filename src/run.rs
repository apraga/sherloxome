//! Generic run types shared across datasets.
//!
//! Output VCF filenames must adhere to this convention for analysis:
//! `SAMPLE-SEQUENCER-CAPTURE-DEPTH{-SILICO}.*`
//! Value are free text for maximum compatibility, with the exception of DEPTH that must be
//! an integer followed by the 'x' character.
//! Here is an example
//!
//!    HG002-novaseq-agilent-100x-simuscop.haplotypecaller.vcf.gz
//!
//! Where SAMPLE is HG002, SEQUENCER is novaseq, DEPTH is 100. It's in silico data as SILICO is simuscop (one of the tools supported, see `crate::silico`). The rest of the filename is discarded.
//! A set of values is defined in `crate::baid2020` for stronger type-checking.
use regex::Regex;
use serde::Deserialize;
use std::path::PathBuf;
use std::sync::OnceLock;

static RUN_RE: OnceLock<Regex> = OnceLock::new();

/// Regexp defining our syntax
fn run_re() -> &'static Regex {
    RUN_RE.get_or_init(|| {
        Regex::new(
            r"^(?P<sample>[^-]+)-(?P<sequencer>[^-]+)-(?P<capture>[^-]+)-(?P<depth>[^-.x]+)x(?:-(?P<silico>[^.-]+))?",
        )
        .unwrap()
    })
}

/// Infer a [`Run`] from a VCF filename following the
/// `SAMPLE-CAPTURE-DEPTH-SEQUENCER{-SILICO}.*.vcf` convention.
pub fn run_from_filename(fname: &PathBuf) -> Option<Run> {
    let name = fname.file_name()?.to_string_lossy();
    let caps = run_re().captures(name.as_ref())?;
    Some(Run {
        sample: caps["sample"].parse().ok()?,
        capture: caps["capture"].parse().ok()?,
        depth: caps["depth"].parse().ok()?,
        sequencer: caps["sequencer"].parse().ok()?,
        silico: caps.name("silico").and_then(|m| m.as_str().parse().ok()),
    })
}

/// Format a [`Run`] as `SAMPLE-SEQUENCER-CAPTURE-DEPTH{-SILICO}`.
pub fn run_to_string(run: &Run) -> String {
    let base = format!(
        "{}-{}-{}-{}x",
        run.sample, run.sequencer, run.capture, run.depth
    );
    match &run.silico {
        None => base,
        Some(s) => format!("{}-{}", base, s),
    }
}

/// A run is defined by a sample (string), [Sequencer], capture [Capture], [Depth] and [Silico] status.
#[derive(Clone, Deserialize, Debug, Hash, Eq, PartialEq)]
#[serde(rename_all = "lowercase")]
pub struct Run {
    pub sample: String,
    pub sequencer: String,
    pub capture: String,
    pub depth: u32,
    pub silico: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::baid2020::{Capture, Sequencer, to_run};
    use crate::giab::Patient;
    use std::path::PathBuf;
    fn make_run(silico: Option<String>) -> Run {
        return to_run(
            Patient::HG002,
            Sequencer::Novaseq,
            Capture::Agilent,
            100,
            silico,
        );
    }

    #[test]
    fn run_to_string_real() {
        assert_eq!(run_to_string(&make_run(None)), "HG002-novaseq-agilent-100x");
    }

    #[test]
    fn run_to_string_silico() {
        assert_eq!(
            run_to_string(&make_run(Some(String::from("simuscop")))),
            "HG002-novaseq-agilent-100x-simuscop"
        );
        assert_eq!(
            run_to_string(&make_run(Some(String::from("varben")))),
            "HG002-novaseq-agilent-100x-varben"
        );
    }

    #[test]
    fn run_from_filename_real() {
        let path = PathBuf::from("HG002-novaseq-agilent-100x.vcf.gz");
        assert_eq!(run_from_filename(&path), Some(make_run(None)));
    }

    #[test]
    fn run_from_filename_silico() {
        let path = PathBuf::from("HG002-novaseq-agilent-100x-simuscop.vcf.gz");
        assert_eq!(
            run_from_filename(&path),
            Some(make_run(Some(String::from("simuscop"))))
        );
    }

    #[test]
    fn run_from_filename_invalid() {
        let bad = [
            "nomatch.vcf",
            "HG002novaseq-bad-100x.vcf",
            "HG002-novaseqbad-100x.vcf",
            "HG002-novaseqbad-100.vcf",
            "HG002-novaseq-bad100x.vcf",
        ];
        for b in bad {
            assert_eq!(run_from_filename(&PathBuf::from(b)), None);
        }
    }

    #[test]
    fn roundtrip() {
        for silico in [
            None,
            Some(String::from("simuscop")),
            Some(String::from("varben")),
        ] {
            let run = make_run(silico);
            let s = run_to_string(&run);
            let path = PathBuf::from(format!("{s}.vcf.gz"));
            assert_eq!(run_from_filename(&path), Some(run));
        }
    }
}
