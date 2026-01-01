## Requirements

Software (TODO flakes)
- hap.py
- vcfeval

Files
- a genome (do a symlink ?)
- sdf from genome
```bash
rtg format /Work/Groups/bisonex/dgenomes/genome-human/GCA_000001405.15_GRCh38_full_analysis_set.fna  -o  sdf
```
-  Download truth files with `make truth`

## Run

Select the combination patient/capture kit/sequecer from the full samplesheet. A run can accomodate several fastq but only for a single capture kit.
For example, select all HG001 data for Agilent and run sarek with

```bash
head -n 1 samplesheet-full.csv  > samplesheet.csv
grep agilent samplesheet-full.csv | grep 'HG001' >> samplesheet.csv
make run-agilent
```

## Capture kits

Baid 2020 : "Agilent v7, IDT-xGen, and Nextera"
