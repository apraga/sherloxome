# Get raw data

Run `sherloxome setup` to download reference data and/or generate in silico controls:

```bash
sherloxome setup -c config.toml
```

This command writes `samplesheet.csv` for the variant calling pipeline.

---

## Real patient data (GIAB)

Add a `[real]` section to `config.toml` to download GIAB (Genome In A Bottle) benchmark data:

```toml
[real]
patients   = ["HG001", "HG002", "HG003"]
sequencers = ["hiseq4000", "novaseq"]
captures   = ["agilent", "idt", "truseq"]
depths     = ["50x", "75x", "100x"]
```

### Available patients

| ID | Sample | Population |
|----|--------|-----------|
| HG001 | NA12878 | CEU (CEPH) |
| HG002 | NA24385 | Ashkenazi son |
| HG003 | NA24149 | Ashkenazi father |
| HG004 | NA24143 | Ashkenazi mother |
| HG005 | NA24631 | Chinese son |
| HG006 | NA24694 | Chinese father |
| HG007 | NA24695 | Chinese mother |

### Available combinations

Not all (patient × sequencer × capture × depth) combinations exist in the BAID2020 dataset. Sherloxome filters your request against the known-available set and prints a summary:

```
You asked for 126 runs (7 patients x 3 kits x 3 depths x 2 sequencers)
Only 84 are available
```

The available combinations are:

- **HiSeq4000 at 50x**: all 7 patients × all 3 kits (Agilent, IDT, TruSeq)
- **NovaSeq**: all 7 patients × Agilent 50x/75x/100x + IDT 50x/75x/100x + TruSeq 50x/75x

### Downloaded files

For each patient, files are saved to `data/ref/`:

- `HG00X_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` — truth VCF
- `HG00X_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi` — tabix index
- `HG00X_GRCh38_1_22_v4.2.1_benchmark[_noinconsistent].bed` — high-confidence BED

> HG002–HG004 use a `_noinconsistent` BED that excludes regions with inconsistent calls across family members.

FASTQ files are **not** downloaded. For real patients, `samplesheet.csv` contains Google Cloud Storage URLs; Nextflow downloads them during the pipeline run.

---

## Synthetic in silico controls

Add a `[silico]` section to inject ClinVar pathogenic variants into a real BAM:

```toml
[silico]
capture     = "agilent"
bam         = "https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam"
clinvar     = "data/exp_raw/clinvar.vcf.gz"   # optional; downloaded from NCBI if absent
nb_variants = 1000
```

### Steps performed

1. The BAM is downloaded (if `bam` is an URL) or verified locally
2. Hard-clipped reads are removed with `samtools` + `awk`
3. `nb_variants` ClinVar pathogenic SNVs are sampled from within the capture BED, enforcing ≥50 bp spacing between selected variants
4. `muteditor` inserts the variants into the BAM with random allele fractions (0.4–0.6)
5. The edited BAM is converted to paired FASTQ via `samtools fastq`
6. A VCF of successfully inserted variants is written alongside the FASTQ

### Variant selection criteria

A ClinVar variant is eligible if it:

- Falls within the capture kit BED intervals
- Has `CLNSIG` of `Pathogenic`, `Likely_pathogenic`, or `Uncertain_significance`
- Is an SNV (single-nucleotide variant, REF and ALT both length 1)
- Is ≥50 bp from the nearest already-selected variant on the same chromosome

### Output files

By default, output is written to `data/exp_raw/`. Override with `outdir`:

```toml
[silico]
outdir = "silico"
```

Key outputs:

| File | Description |
|------|-------------|
| `clinvar_{capture}.vcf.gz` | Sampled ClinVar variants (truth VCF) |
| `clinvar_{capture}.mut` | Mutation file for muteditor |
| `varben/edit.sorted.bam` | BAM with variants injected |
| `{sample}_success.vcf.gz` | VCF of successfully inserted variants |
| `{sample}_1.fq.gz` / `{sample}_2.fq.gz` | Paired FASTQ ready for the pipeline |

### Reference genome

Set the reference FASTA path in `config.toml`:

```toml
fasta = "data/ref/GCA_000001405.15_GRCh38_full_analysis_set.fna"
```

If the path does not exist, sherloxome downloads the GRCh38 full analysis set from NCBI (~3 GB compressed) and decompresses it. The BWA index (~10 GB) is also downloaded automatically.
