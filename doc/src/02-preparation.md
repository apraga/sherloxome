# Prepare FASTQ

## Quick start 
To download reference data and/or generate in silico controls, edit `config.toml` to select which step you want.
Then run the CLI with `setup`. We suggest using nix:

```bash
nix shell .#default --command cargo run --release -- setup
```

This command writes one `samplesheet-{capture}.csv` per capture kit for the variant calling pipeline. Splitting by capture kit is necessary because a single sarek run only accepts one `--intervals` BED, so rows using different capture kits cannot share a samplesheet.

## Detailed example

Here's an example combining 1 real patient (capture `agilent`) and 2 insilico configurations (capture `agilent-col6a1`). Two files are generated:

`samplesheet-agilent.csv`:
```csv
patient,sample,lane,fastq_1,fastq_2
HG002,HG002_hiseq4000_agilent_50x,1,https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes_agilent.50x.R1.fastq.gz,https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes_agilent.50x.R2.fastq.gz
```

`samplesheet-agilent-col6a1.csv`:
```csv
patient,sample,lane,fastq_1,fastq_2
silico-varben,HG002_hiseq4000_agilent-col6a1_50x_varben_varben,1,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_varben_1.fq.gz,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_varben_2.fq.gz
silico-simuscop,nopatient_hiseq4000_agilent-col6a1_50x_simuscop_simuscop,1,data/exp_raw/nopatient_hiseq4000_agilent-col6a1_50x_simuscop_1.fq.gz,data/exp_raw/nopatient_hiseq4000_agilent-col6a1_50x_simuscop_2.fq.gz
```

And the correspondig `config.toml`:

```toml
 # URL to fasta, or link to local version, otherwise download it from NCBI
fasta = "data/ref/GCA_000001405.15_GRCh38_full_analysis_set.fna"

# Configuration for raw GIAB patients data
[real]
patients = [HG002] 
depths = [50] 
captures = ["agilent"] 
sequencers = ["hiseq4000"] 

[silico]
# VCF containing clinvar variants. 
clinvar = "data/exp_raw/clinvar.vcf.gz"
# Number of random clinvar variants to insert into the BAM file. Default is 1000
nb_variants = 1000
# Capture kit name (the BED must be defined in [capture] below)
capture = "agilent"
# A BAM file is always required for varben as it will will modify it
bam_file = "data/exp_raw/HG002_hiseq4000_agilent_50x.bam"

# This section enables simuscop FASTQ generation
[silico.simuscop]
# Pre-built seqToProfile profile directory. 
profile = "data/ref/hiseq4000_agilent_50x.profile"
# Sequencing coverage
coverage = 50

# This section enables varben BAM editing
[silico.varben]
mindepth = 50

# Define here the name of alls captures and the bed file
[capture]
# Those values are mandatory for GIAB data (Baid 2020)
idt = "data/capture/idt_capture.grch38.bed"
truseq = "data/capture/truseq-dna-exome-targeted-regions-manifest-v1-2-lifted-grch38.bed"
agilent = "data/capture/agilent.targets.grch38.bed"
# For testing
agilent-col6a1 = "data/capture/agilent-col6a1.targets.grch38.bed"
```

See also the [filenaming scheme](050-filenaming.md).

