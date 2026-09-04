# FAQ & Troubleshooting

## hap.py fails with a vcfeval error about awk, hostname not found, or too many threads

vcfeval has an undocumented thread limit. Sherloxome forces `--threads 1` for each hap.py call. If you are calling hap.py manually, explicitly set `--threads 1`.

In a SLURM job the `SLURM_CPUS_PER_TASK` variable can cause hap.py to request more threads than vcfeval allows. Using sherloxome's `analyze` command avoids this.

## Only some of my requested runs appear in the generated samplesheets

Not all (patient × sequencer × capture × depth) combinations are available in the BAID2020 dataset. Sherloxome prints a summary of requested vs. available runs during `setup`. See [available combinations](04-misc.md) for the full list.

## The BAM download stalls or is very slow

Large BAM files can be several GB. The HTTP client has a 30-second connection timeout but no total transfer timeout, so the download will not abort mid-transfer. Use `RUST_LOG=info` to see progress.

## ClinVar is missing

If `clinvar` is not set in `config.toml`, sherloxome downloads it from NCBI to `data/exp_raw/clinvar.vcf.gz`. To reuse a local copy:

```toml
[silico]
clinvar = "/path/to/clinvar.vcf.gz"
```

## The reference FASTA is missing

If the `fasta` path does not exist, sherloxome downloads the GRCh38 full analysis set from NCBI (~3 GB compressed) and decompresses it to `data/ref/`. The BWA index (~10 GB) is also downloaded automatically on first use.

## muteditor is very slow

`muteditor` realigns reads around each injected variant using BWA. For 1000 variants this can take several hours. Run on a cluster — see the [SLURM example](04-misc.md).

## `sherloxome analyze` skips all VCFs

VCF filenames must contain the patient ID, sequencer, capture kit, and depth as they appear in their string representations (e.g. `HG002`, `novaseq`, `agilent`, `75x`). Sarek names output VCFs after the `sample` column of the samplesheet, which sherloxome formats correctly. If you renamed files manually, check that all four fields are present in each filename.

## The plot window does not open

`sherloxome plot` calls the system browser via the Vega-Lite library. Ensure a graphical session is available. On a remote server, copy `merged.csv` to your local machine and run `plot` there.
