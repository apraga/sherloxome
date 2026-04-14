# Overview

1. Download reference data (VCF + BED) : `sherloxome download`
2. Generate or download fastq
  - if you use sequencing data that are avaiable, online: TODO
  - for in silico data (or real data that is modified ) : TODO
3. Generate samplesheet (nf-core/sarek only)
4. Run the pipeline
5. Analyse the result `sherloxome analyze`

## Requirements

Software (TODO flakes)
- hap.py
- vcfeval

Rust (rustup for now)

## Download real patients data (GIAB)

Generate a samplesheet `samplesheet.csv` with

```bash
cargo run --release -c config.toml samplesheet
```

We have to run the pipeline for each kit. So The following configuration file will generate all combinations
for agilent
```toml
[real]
patients = ["HG001", "HG002", "HG003", "HG004", "HG005", "HG006", "HG007"]
depths = ["50x", "75x", "100x"]
kits = ["agilent"]
sequencers = ["hiseq4000", "novaseq"]

[silico]
patients = ["HG001"]
```


## Generate controls

Setup your config.toml to use a BAM file
```toml
[silico]
capture = "agilent"
# Either an URL or a local link to insert variants in a BAM file
bam = "https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam"
# Required for generating a bam file
fasta = "/Work/Groups/bisonex/dgenomes/genome-human/GCA_000001405.15_GRCh38_full_analysis_set.fna"
```

Then run
```bash
cargo run --release controls
```
This requires bwa, samtools and varben to be installed beforehand. `varben` is rather long so it should be started from a slurm job directly.
TODO : what this does

## Run

Assuming you want to run all agilent runs, and that you use repro-sarek, go into repro-sarek directory, then

```bash
nix develop
```
Now go into sherloxome
```bash
   nextflow run nf-core/sarek --input samplesheet.csv -r 3.5.1 --outdir data/exp_raw/giab  -c ../repro-sarek/nextflow.config  -c ../repro-sarek/conf/slurm.config --tools haplotypecaller --skip_tools haplotypecaller_filter --wes --intervals data/capture/Agilent_SureSelect_All_Exons_v7_hg38_Regions.bed
````

## Configuration

Configuration is done in a `config.toml` defining the data to be analyzed. For "real" patients, define which GIAB patient it is, which capture kitand which sequencer is used.
