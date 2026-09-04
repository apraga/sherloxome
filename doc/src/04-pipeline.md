# Run the variant calling pipeline

After `sherloxome setup` writes one `samplesheet-{capture}.csv` per capture kit, run a variant calling pipeline on the data.

## nf-core/sarek

Sherloxome generates samplesheets compatible with [nf-core/sarek](https://nf-co.re/sarek).

### Samplesheet format

`setup` groups rows by capture kit and writes a separate file for each, e.g. `samplesheet-agilent.csv`:

```csv
patient,sample,lane,fastq_1,fastq_2
HG001,HG001_hiseq4000_agilent_50x,1,https://...R1.fastq.gz,https://...R2.fastq.gz
HG002,HG002_novaseq_agilent_75x,1,https://...R1.fastq.gz,https://...R2.fastq.gz
silico-varben,HG002_nohardclip,1,/path/to/HG002_nohardclip_1.fq.gz,/path/to/HG002_nohardclip_2.fq.gz
```

- Real patient rows contain GCS FASTQ URLs; Nextflow downloads them at runtime.
- In silico rows contain local FASTQ paths produced by `sherloxome setup`.
- Splitting by capture kit is required because a single sarek run only accepts one `--intervals` BED, so samples using different capture kits cannot share a run.

### Example sarek command

Run sarek once per capture kit, using the matching samplesheet and `--intervals`:

```bash
nextflow run nf-core/sarek \
    --input samplesheet-agilent.csv \
    -r 3.5.1 \
    --outdir data/exp_raw/giab \
    --tools haplotypecaller \
    --skip_tools haplotypecaller_filter \
    --wes \
    --intervals data/capture/agilent.targets.grch38.bed
```

### Capture kit BED files

Configure paths in the `[capture]` section of `config.toml`:

```toml
[capture]
agilent = "data/capture/agilent.targets.grch38.bed"
idt     = "data/capture/idt_capture.grch38.bed"
truseq  = "data/capture/truseq-dna-exome.bed"
```

| Key | BED source |
|-----|-----------|
| `agilent` | Agilent SureSelect All Exons v7 (hg38) |
| `idt` | IDT xGen Exome Hyb Panel v2 (hg38) |
| `truseq` | Illumina TruSeq DNA Exome (lifted to hg38 via UCSC) |

## Output VCF naming

The `analyze` step infers run metadata (patient, sequencer, capture, depth) from VCF filenames. Sarek names output VCFs using the `sample` column of the samplesheet, so the naming is already correct as long as you use the samplesheet written by `sherloxome setup`.

Expected filename pattern: `{patient}_{sequencer}_{capture}_{depth}x.vcf.gz`

Example: `HG002_novaseq_agilent_75x.vcf.gz`

## SLURM example

```bash
#!/bin/bash
#SBATCH --job-name=sarek
#SBATCH --output=%x.%J.out
#SBATCH --time=24:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

nextflow run nf-core/sarek \
    --input samplesheet-agilent.csv \
    -r 3.5.1 \
    --outdir data/exp_raw/giab \
    --tools haplotypecaller \
    --skip_tools haplotypecaller_filter \
    --wes \
    --intervals data/capture/agilent.targets.grch38.bed
```
