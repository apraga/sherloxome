# Sherloxome

[![CI](https://github.com/apraga/sherloxome/actions/workflows/ci.yml/badge.svg)](https://github.com/apraga/sherloxome/actions/workflows/ci.yml)

A benchmarking tool for validating exome/targeted capture. It serves to prepare testing FASTQ data (real patients or insilico) and compute benchmarking data.

## Quickstart

Install all dependencies with Nix:
```bash
nix profile add .
```
Generate FASQT for testing:
```bash
cargo install
sherloxome setup
```
Run your pipeline. For germline analysis, we recommend [nf-core/sarek](https://nf-co.re/sarek/) with [our reproducible setup](https://github.com/apraga/reproducible-sarek-germline/). After starting a shell with `nix develop`:
```bash
/nextflow run  nf-core/sarek --input samplesheet.csv -r 3.7.1
| --outdir data/exp_raw/giab --tools haplotypecaller --skip_tools haplotypecaller_filter --wes --intervals data/capture/Agilent_SureSelect_All_Exons_v7_hg38_Regions.bed -profile apptainer
```

Analyse data
```bash
sherloxome benchmark
```
