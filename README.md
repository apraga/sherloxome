## Requirements

Software (TODO flakes)
- hap.py
- vcfeval

Rust (rustup for now)

## Analyse real patients

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

### Run
Assuming you want to run all agilent runs, and that you use repro-sarek, go into repro-sarek directory, then

```bash
nix develop
```
Now go into sherloxome
```bash
   nextflow run nf-core/sarek --input samplesheet.csv -r 3.5.1 --outdir data/exp_raw/giab  -c ../repro-sarek/nextflow.config  -c ../repro-sarek/conf/slurm.config --tools haplotypecaller --skip_tools haplotypecaller_filter --wes --intervals data/capture/Agilent_SureSelect_All_Exons_v7_hg38_Regions.bed
````
