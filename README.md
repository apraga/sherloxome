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
```

## Run

With sarek

nextflow run nf-core/sarek --input samplesheet-full.csv  -r 3.5.1 --outdir hg002 --tools haplotypecaller
