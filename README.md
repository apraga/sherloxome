# Sherloxome

A benchmarking tool for validating exome/targeted capture. It serves to prepare testing FASTQ data (real patients or insilico) and compute benchmarking data.

## Quickstart

Generate FASQT for testing:
```bash
cargo install
sherloxome setup
```
Run your pipeline. For `nf-core/sarek`, we recommend [our setup for germline analysis](https://github.com/apraga/reproducible-sarek-germline/):
```bash
/nextflow run  nf-core/sarek --input samplesheet.csv -r 3.7.1
| --outdir data/exp_raw/giab --tools | haplotypecaller --skip_tools haplotypecaller_filter --wes --intervals | data/capture/Agilent_SureSelect_All_Exons_v7_hg38_Regions.bed -profile | apptainer
```

Analyse data
```bash
sherloxome benchmark
```
