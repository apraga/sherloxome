# Sherloxome

![Real patients](https://github.com/apraga/sherloxome/actions/workflows/giab-col6a1.yml/badge.svg)
![Silico (varben)](https://github.com/apraga/sherloxome/actions/workflows/varben-col6a1.yml/badge.svg)

A benchmarking tool for validating exome/targeted capture. It serves to prepare testing FASTQ data (real patients or insilico) and compute benchmarking data.

## Quickstart

Install all dependencies with [Nix](https://nixos.org/download/#download-nix) and our CLI tool with [Rust](https://rustup.rs/):
```bash
nix profile add .
cargo install
```
Generate FASTQ for testing:
```bash
sherloxome setup
```
Run your pipeline. `setup` writes one `samplesheet-{capture}.csv` per capture kit, since sarek only takes a single `--intervals` BED per run. For germline analysis, we recommend [nf-core/sarek](https://nf-co.re/sarek/) with [our reproducible setup](https://github.com/apraga/reproducible-sarek-germline/). After starting a shell with `nix develop`:
```bash
nextflow run  nf-core/sarek --input samplesheet-agilent.csv -r 3.7.1 --outdir data/exp_raw/giab --tools haplotypecaller --skip_tools haplotypecaller_filter --wes --intervals data/capture/Agilent_SureSelect_All_Exons_v7_hg38_Regions.bed -profile apptainer
```

Analyse data:
```bash
sherloxome benchmark
```
