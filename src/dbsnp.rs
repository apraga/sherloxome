//! # Subset of dbSNP for simuscop
//!
//! Select dbSNP variants that are in a capture kit, with AF >= 1% and not in our clinvar sampling.
//! Used by [`crate::simuscop`]).
//!
//! dbSNP uses RefSeq accessions for chromosome names (e.g. `NC_000021.9` for chr21), so a
//! mapping file is required. See the [dbSNP setup
//! guide](https://apraga.github.io/sherloxome/022-dbsnp.html#configuration).
//!
//! dbSNP is never downloaded completely, only data in the capture kit, and the rest (chromosome
//! renaming, AF filtering, excluding clinvar) is delegated to `bcftools` rather than parsed here.
use std::error::Error;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::silico::index_vcf;

/// Latest NCBI dbSNP GRCh38
// const DBSNP_URL: &str = "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz";

/// Run a shell pipeline, e.g. "cmd1 | cmd2 > out"
fn run(cmd: &str) -> Result<(), Box<dyn Error>> {
    let status = Command::new("sh").args(["-c", cmd]).status()?;
    if !status.success() {
        return Err(format!("Command failed ({status}): {cmd}").into());
    }
    Ok(())
}

/// Query dbSNP restricted to the capture kit, keep common SNVs absent from `clinvar_vcf`
/// dbSNP has already been processed to filter only common SNV in the capture kit.
/// To use it with your own capture kit, you need the following command:
/// ```bash
///  # 1. capture bed: chr -> RefSeq (for `bcftools view -R`)
///  awk -F'\t' 'NR==FNR{gsub(/\r/,""); m[$1]=$2; next} $1 in m{print m[$1]"\t"$2"\t"$3}' \
///      "$MAPPING" "$BED" > ${OUT}.regions.bed
///  # 2. build the reverse map: RefSeq -> chr (for `bcftools annotate --rename-chrs`)
///  awk -F'\t' '{gsub(/\r/,""); print $2"\t"$1}' "$MAPPING" > ${OUT}.rename.txt
///  # 3. filter dbSNP to the capture kit, common variants, SNVs only
///  bcftools view -v snps -i 'COMMON=1' -R ${OUT}.regions.bed "$DBSNP" -Oz -o ${OUT}.refseq.vcf.gz
///  # 4. convert result back to chr notation
///  bcftools annotate --rename-chrs ${OUT}.rename.txt ${OUT}.refseq.vcf.gz -Oz -o ${OUT}.vcf.gz
///  # 5. index
///  bcftools index -t ${OUT}.vcf.gz
/// ```bash
pub fn sample_dbsnp(
    capture: &str,
    clinvar_vcf: &Path,
    outdir: &Path,
) -> Result<PathBuf, Box<dyn Error>> {
    let vcf_out = PathBuf::from(format!("data/exp_raw/dbsnp_{capture}.vcf.gz"));
    if vcf_out.exists() {
        log::info!("Skipping dbsnp querying as {:?} already exists", vcf_out);
        Ok(vcf_out)
    } else {
        let common = PathBuf::from(format!("data/exp_raw/dbsnp_{capture}_common.vcf.gz"));
        if !common.exists() {
            let msg = format!("{} should exist. The version for agilent, IDT and truseq have been shipped in the code in data/exp_raw.
            If you are using another capture kit, see the documentation to filter dbSNP data", common.display());
            log::error!("{msg}");
            Err(msg.into())
        } else {
            std::fs::create_dir_all(outdir)?;
            run(&format!(
                "bcftools isec -C -w1 -Oz -o {} {} {}",
                vcf_out.display(),
                common.display(),
                clinvar_vcf.display(),
            ))?;
            index_vcf(&vcf_out)?;

            Ok(vcf_out)
        }
    }
}

/// Write simuscop's SNP background format (see SimuSCoP_User_Guide.pdf §4.2.2):
/// `rsID  chrom  pos  REF/ALT  strand  REF`. dbSNP GRCh38 coordinates are always on the
/// forward (+) strand.
pub fn write_snp_input(dbsnp_vcf: &Path, out: &Path) -> Result<(), Box<dyn Error>> {
    run(&format!(
        "bcftools query -f '%ID\\t%CHROM\\t%POS\\t%REF/%ALT\\t+\\t%REF\\n' {} > {}",
        dbsnp_vcf.display(),
        out.display(),
    ))
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use std::fs;

//     fn write_mapping(dir: &Path) -> PathBuf {
//         let path = dir.join("mapping.txt");
//         // CRLF on purpose: mirrors the upstream NCBI assembly report, which the doc-documented
//         // curl+awk command carries a stray \r from.
//         fs::write(&path, "chr21\r\tNC_000021.9\r\nchr22\tNC_000022.11\n").unwrap();
//         path
//     }

//     /// Bgzip + tabix-index a plain-text VCF, as `bcftools`/`tabix` expect.
//     fn index(dir: &Path, name: &str, content: &str) -> PathBuf {
//         let plain = dir.join(format!("{name}.vcf"));
//         fs::write(&plain, content).unwrap();
//         let gz = dir.join(format!("{name}.vcf.gz"));
//         let status = Command::new("bcftools")
//             .args(["view", "-Oz", "-o"])
//             .arg(&gz)
//             .arg(&plain)
//             .status()
//             .unwrap();
//         assert!(status.success());
//         index_vcf(&gz).unwrap();
//         gz
//     }

//     #[test]
//     fn prepare_regions_and_rename_map_translates_and_skips_unknown_chroms() {
//         let dir = std::env::temp_dir().join("dbsnp_prepare_test");
//         fs::create_dir_all(&dir).unwrap();

//         let bed = dir.join("capture.bed");
//         fs::write(&bed, "chr21\t100\t200\nchrUn_foo\t0\t10\nchr22\t300\t400\n").unwrap();

//         let regions = dir.join("regions.bed");
//         let rename = dir.join("rename.txt");
//         prepare_regions_and_rename_map(&write_mapping(&dir), &bed, &regions, &rename).unwrap();

//         assert_eq!(
//             fs::read_to_string(&regions).unwrap(),
//             "NC_000021.9\t100\t200\nNC_000022.11\t300\t400\n"
//         );
//         assert_eq!(
//             fs::read_to_string(&rename).unwrap(),
//             "NC_000021.9\tchr21\nNC_000022.11\tchr22\n"
//         );
//     }

//     #[test]
//     fn write_snp_input_matches_simuscop_format() {
//         let dir = std::env::temp_dir().join("dbsnp_snp_input_test");
//         fs::create_dir_all(&dir).unwrap();

//         let vcf = index(
//             &dir,
//             "variants",
//             "##fileformat=VCFv4.2\n\
//              ##contig=<ID=chr21>\n\
//              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
//              chr21\t100\trs123\tA\tG\t.\t.\t.\n",
//         );

//         let out = dir.join("out.snp");
//         write_snp_input(&vcf, &out).unwrap();

//         assert_eq!(
//             fs::read_to_string(&out).unwrap(),
//             "rs123\tchr21\t100\tA/G\t+\tA\n"
//         );
//     }

//     /// End-to-end pipeline against local fixtures (no network): capture kit restriction,
//     /// chromosome renaming, multiallelic split, COMMON/SNV filtering, and ClinVar exclusion.
//     #[test]
//     fn sample_dbsnp_filters_and_renames_locally() {
//         let dir = std::env::temp_dir().join("dbsnp_sample_test");
//         let _ = fs::remove_dir_all(&dir);
//         fs::create_dir_all(&dir).unwrap();

//         let header = "##fileformat=VCFv4.2\n\
//              ##contig=<ID=NC_000021.9>\n\
//              ##INFO=<ID=COMMON,Number=0,Type=Flag,Description=\"common\">\n\
//              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
//         let dbsnp_vcf = index(
//             &dir,
//             "dbsnp",
//             &format!(
//                 "{header}\
//                  NC_000021.9\t100\trs1\tA\tG,T\t.\t.\tCOMMON\n\
//                  NC_000021.9\t150\trs2\tC\tT\t.\t.\t.\n\
//                  NC_000021.9\t200\trs3\tAT\tA\t.\t.\tCOMMON\n\
//                  NC_000021.9\t300\trs4\tG\tC\t.\t.\tCOMMON\n\
//                  NC_000021.9\t999\trs5\tA\tG\t.\t.\tCOMMON\n"
//             ),
//         );
//         let clinvar_vcf = index(
//             &dir,
//             "clinvar",
//             "##fileformat=VCFv4.2\n##contig=<ID=chr21>\n\
//              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
//              chr21\t300\t.\tG\tC\t.\t.\t.\n",
//         );

//         let mapping = write_mapping(&dir);
//         let bed = dir.join("capture.bed");
//         fs::write(&bed, "chr21\t50\t250\n").unwrap(); // covers rs1..rs3, excludes rs4/rs5

//         let config = DbsnpConfig {
//             vcf: Some(dbsnp_vcf.to_string_lossy().to_string()),
//             mapping: Some(mapping),
//         };

//         let out = sample_dbsnp(&config, &bed, "test", &clinvar_vcf, &dir).unwrap();
//         let output = Command::new("bcftools")
//             .args(["query", "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\n"])
//             .arg(&out)
//             .output()
//             .unwrap();
//         let content = String::from_utf8(output.stdout).unwrap();

//         // rs1 (COMMON, multiallelic) is split and kept; rs2 (not COMMON) is dropped; rs3
//         // (COMMON, indel) is dropped by the SNV filter; rs4 (COMMON, at the ClinVar position)
//         // is dropped by isec; rs5 is outside the capture kit.
//         assert_eq!(
//             content,
//             "chr21\t100\trs1\tA\tG\n\
//              chr21\t100\trs1\tA\tT\n"
//         );
//     }
// }
