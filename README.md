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
<!-- Files -->

<!-- - a genome (do a symlink ?) -->
<!-- - sdf from genome -->
<!-- ```bash -->
<!-- rtg format /Work/Groups/bisonex/dgenomes/genome-human/GCA_000001405.15_GRCh38_full_analysis_set.fna  -o  sdf -->
<!-- ``` -->
<!-- -  Download truth files with `make truth` -->

<!-- ## Run -->

<!-- Select the combination patient/capture kit/sequecer from the full samplesheet. A run can accomodate several fastq but only for a single capture kit. -->
<!-- For example, select all HG001 data for Agilent and run sarek with -->

<!-- ```bash -->
<!-- head -n 1 samplesheet-full.csv  > samplesheet.csv -->
<!-- grep agilent samplesheet-full.csv | grep 'HG001' >> samplesheet.csv -->
<!-- make run-agilent -->
<!-- ``` -->

<!-- ## Capture kits -->

<!-- Baid 2020 : "Agilent v7, IDT-xGen, and Nextera" -->

<!-- ## Analysis -->

<!-- `make compare` will run  `compare.sh`, which will run hap.py over all vcf genereted by sarek in `giab` (path is harcoded). By default, parallel uses 4 cores, 1 for each hap.py. -->
<!-- All output files will be is the `analysis` folder. -->

<!-- `make merge` concanate and reformat the summarized output of all happy comparisons into `analysis/all-giab.csv` -->
