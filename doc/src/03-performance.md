# Evaluate performance

## Run hap.py analysis

After variant calling, benchmark the output VCFs against GIAB truth sets:

```bash
sherloxome analyze \
    -i data/exp_raw/giab \
    -o data/analysis \
    -c config.toml
```

| Flag | Description |
|------|-------------|
| `-i` / `--input` | Directory to search recursively for `**/*.vcf.gz` |
| `-o` / `--output` | Output directory for hap.py summaries and `merged.csv` |
| `-c` / `--config` | Config file (needed for capture BED paths; default `config.toml`) |

## What `analyze` does

For each VCF whose filename contains recognisable run metadata (patient, sequencer, capture, depth):

1. Locates the GIAB truth VCF and high-confidence BED in `data/ref/`
2. Locates the capture kit BED from `[capture]` in `config.toml`
3. Generates an RTG SDF from the reference FASTA on first run (`rtg format`)
4. Runs `hap.py` with the `vcfeval` engine (single-threaded per run; VCFs are processed in parallel)
5. Writes per-run hap.py summaries to the output directory
6. Merges all summaries into `merged.csv`, adding `patient`, `capture`, `sequencer`, `depth` columns

## Output files

```
data/analysis/
├── HG001-agilent-hiseq4000-50x.summary.csv
├── HG002-agilent-novaseq-75x.summary.csv
├── ...
└── merged.csv
```

`merged.csv` is the input to the `plot` step.

## Prerequisites

- GIAB truth VCFs in `data/ref/` (written by `sherloxome setup`)
- The GRCh38 reference FASTA specified in `config.toml`
- `hap.py` and `rtg` in `PATH`

## Restart safety

Any run for which a `.summary.csv` already exists in the output directory is skipped. You can safely re-run `analyze` after adding new VCFs without reprocessing existing results.

## Visualising results

Once `merged.csv` exists, open an interactive boxplot:

```bash
sherloxome plot -i data/analysis/merged.csv
```

The chart shows F1-score distributions broken down by:

- Patient (HG001–HG007)
- Capture kit (Agilent, IDT, TruSeq)
- Sequencer (HiSeq4000, NovaSeq)
- Depth (50x, 75x, 100x)

Only `PASS` variants are included in the plot.

## SLURM example

```bash
#!/bin/bash
#SBATCH --job-name=analyze
#SBATCH --output=%x.%J.out
#SBATCH --time=4:00:00
#SBATCH --partition=smp
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G

sherloxome analyze \
    -i data/exp_raw/giab \
    -o data/analysis \
    -c config.toml
```
