# Miscellaneous

## Logging

Control verbosity with the `RUST_LOG` environment variable:

```bash
RUST_LOG=info  sherloxome setup    # progress and download messages
RUST_LOG=debug sherloxome setup    # detailed trace: file paths, skipped steps, etc.
```

## Restart safety

All long-running steps check whether their output already exists before running:

| Step | Skipped if |
|------|-----------|
| File download | Output file already exists |
| ClinVar sampling | Sampled VCF and mutation file both exist |
| `muteditor` BAM editing | `varben/{sample}/edit.sorted.bam` exists |
| FASTQ generation | Both `_1.fq.gz` and `_2.fq.gz` exist |
| hap.py analysis | `.summary.csv` exists for that run |

This means you can safely interrupt and restart any step.

## Supported capture kits

| Config key | Full name |
|-----------|-----------|
| `agilent` | Agilent SureSelect All Exons v7 (hg38) |
| `idt` | IDT xGen Exome Hyb Panel v2 (hg38) |
| `truseq` | Illumina TruSeq DNA Exome (lifted to hg38 via UCSC) |

## Available run combinations

Not all (patient × sequencer × capture × depth) combinations exist in the BAID2020 dataset. The table below summarises what is available.

**HiSeq4000 at 50x** — all 7 patients × all 3 kits:

| Patient | Agilent 50x | IDT 50x | TruSeq 50x |
|---------|:-----------:|:-------:|:----------:|
| HG001–HG007 | ✓ | ✓ | ✓ |

**NovaSeq** — all 7 patients:

| Capture | 50x | 75x | 100x |
|---------|:---:|:---:|:----:|
| Agilent | ✓ | ✓ | ✓ |
| IDT | ✓ | ✓ | ✓ |
| TruSeq | ✓ | ✓ | ✗ |

## Running `muteditor` on a cluster

`muteditor` performs per-variant BWA realignment and is slow for large variant sets. For `nb_variants = 1000`, expect several hours on a single node. Submit `sherloxome setup` as a SLURM job:

```bash
#!/bin/bash
#SBATCH --job-name=sherloxome-setup
#SBATCH --output=%x.%J.out
#SBATCH --time=8:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G

sherloxome setup -c config.toml
```

## Running the steps separately

`sherloxome setup` runs 4 steps in a row. To
| Step | Runs | Description |
|------|------|-------------|
| `setup prepare` | once | Download reference data, sample clinvar and dbSNP variants |
| `setup simuscop --profile P` | per profile | FASTQ of one simuscop profile. Capture kit and depth come from the filename `SEQUENCER_CAPTURE_DEPTHx.profile` |
| `setup varben --bam B` | per BAM | Insert the variants in one BAM and write its FASTQ. Capture kit comes from the BAM filename |
| `setup samplesheet` | once | Write `samplesheet-{capture}.csv` from every silico FASTQ in the output directory |

To run setup for multiple configuration, for example with a Slurm array:

- `prepare` must run once to sample variants only once
- `simuscop` and `varben` can be run for multiple profiles. You can put them in a Slurm array. They will fail with an error instead of sampling if they are missing.
- `samplesheet` will merge all fastq generated in a single samplesheet

To run multiple simuscop with a slurm array, create a list of profiles

Then run `sbatch` on this script

```
  
#!/bin/bash -l

#SBATCH --job-name="simuscop"
#SBATCH -t 04:00:00
#SBATCH --partition=???
#SBATCH -c 4  ## request 16 cores (MAX is 32)
#SBATCH --mem=12G ## (MAX is 96G)
#SBATCH --array=1-15
sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/profiles.txt
./sherloxome setup simuscop --profile $(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/profiles.txt)
#+end_src
```
