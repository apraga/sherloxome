# Real patient data (GIAB)

Add a `[real]` section to `config.toml` to download GIAB (Genome In A Bottle) benchmark data:

```toml
[real]
patients   = ["HG001", "HG002", "HG003"]
sequencers = ["hiseq4000", "novaseq"]
captures   = ["agilent", "idt", "truseq"]
depths     = ["50x", "75x", "100x"]
```

## Available patients

| ID | Sample | Population |
|----|--------|-----------|
| HG001 | NA12878 | CEU (CEPH) |
| HG002 | NA24385 | Ashkenazi son |
| HG003 | NA24149 | Ashkenazi father |
| HG004 | NA24143 | Ashkenazi mother |
| HG005 | NA24631 | Chinese son |
| HG006 | NA24694 | Chinese father |
| HG007 | NA24695 | Chinese mother |

## Available combinations

Not all (patient × sequencer × capture × depth) combinations exist in the BAID2020 dataset. Sherloxome filters your request against the known-available set and prints a summary:

```
You asked for 126 runs (7 patients x 3 kits x 3 depths x 2 sequencers)
Only 84 are available
```

The available combinations are:

- **HiSeq4000 at 50x**: all 7 patients × all 3 kits (Agilent, IDT, TruSeq)
- **NovaSeq**: all 7 patients × Agilent 50x/75x/100x + IDT 50x/75x/100x + TruSeq 50x/75x

## Downloaded files

For each patient, files are saved to `data/ref/`:

- `HG00X_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` — truth VCF
- `HG00X_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi` — tabix index
- `HG00X_GRCh38_1_22_v4.2.1_benchmark[_noinconsistent].bed` — high-confidence BED

> HG002–HG004 use a `_noinconsistent` BED that excludes regions with inconsistent calls across family members.

FASTQ files are **not** downloaded. For real patients, `samplesheet.csv` contains Google Cloud Storage URLs; Nextflow downloads them during the pipeline run.

