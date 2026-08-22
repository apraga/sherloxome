# Synthetic in silico controls

Silico data consits of a FAST generated either
- from a patient BAM with injected Clinvar pathogenic variants with `varben`
- purely in silico for a sequencer and capture kit based on a model with `simuscop`
Both can be combined. See [Limitations](07-limitations.md)

## Variant selection criteria

A ClinVar variant is eligible if it:

- Falls within the capture kit BED intervals
- Has `CLNSIG` of `Pathogenic`, `Likely_pathogenic`, or `Uncertain_significance`
- Is an SNV (single-nucleotide variant, REF and ALT both length 1)
- Is ≥50 bp from the nearest already-selected variant on the same chromosome


## Common configuration

A `[silico]` section contains setup for both `simuscop` and `varben` for clinvar variants insertion :

```toml
[silico]
 VCF containing clinvar variants. If not set, the VCF will be downloaded from NCBI
clinvar = "data/exp_raw/clinvar_col6a1.vcf.gz"
 Number of random clinvar variants to insert into the BAM file. Default is 1000
nb_variants = 2
 URL to fasta, or link to local version, otherwise download it from NCBI
 fasta =  "https://github.com/nf-core/test-datasets/blob/sarek3/data/genomics/homo_sapiens/genome/chr21/sequence/genome.fasta"
```

and a BAM file for varben is also needd in the common section as it can be used by simuscop or varben, but for different purposes.

```toml
 Capture kit name (the BED must be defined in [capture] below)
capture = "agilent"
 A BAM file is always required for varben as it will will modify it
 Simuscop will define a sequencing profile from it if it does not exist.
 Can be a URL or a local path.
bam_file = "data/exp_raw/HG002.hiseq4000.wes-agilent.50x.dedup.grch38_nohardclip.bam"
```

