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

```toml
 Capture kit name (the BED must be defined in [capture] below)
capture = "agilent"
```

Simuscop needs a pre-built seqToProfile profile, in `[silico.simuscop] profile` — see
[Silico FASTQ with simuscop](0211-simuscop.md). Varben has its own list of BAMs to edit, in
`[silico.varben] bam_files` — see [Varben](0210-varben.md).

