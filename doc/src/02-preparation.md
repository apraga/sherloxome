# Get raw data

Run `sherloxome setup` to download reference data and/or generate in silico controls:

```bash
sherloxome setup -c config.toml
```

This command writes `samplesheet.csv` for the variant calling pipeline. Here's an example combining 1 real patient and 2 insilico configurations

```csv
patient,sample,lane,fastq_1,fastq_2
HG002,HG002_hiseq4000_agilent_50x,1,https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes_agilent.50x.R1.fastq.gz,https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes_agilent.50x.R2.fastq.gz
silico-varben,HG002_hiseq4000_agilent-col6a1_50x_varben_varben,1,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_varben_1.fq.gz,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_varben_2.fq.gz
silico-simuscop,nopatient_hiseq4000_agilent-col6a1_50x_simuscop_simuscop,1,data/exp_raw/nopatient_hiseq4000_agilent-col6a1_50x_simuscop_1.fq.gz,data/exp_raw/nopatient_hiseq4000_agilent-col6a1_50x_simuscop_2.fq.gz
```

See also the [filenaming scheme](050-filenaming.md).
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

Silico data consits of a FAST generated either
- from a patient BAM with injected Clinvar pathogenic variants with `varben`
- purely in silico for a sequencer and capture kit based on a model with `simuscop`
Both can be combined. See [Limitations](07-limitations.md)

### Variant selection criteria

A ClinVar variant is eligible if it:

- Falls within the capture kit BED intervals
- Has `CLNSIG` of `Pathogenic`, `Likely_pathogenic`, or `Uncertain_significance`
- Is an SNV (single-nucleotide variant, REF and ALT both length 1)
- Is ≥50 bp from the nearest already-selected variant on the same chromosome


### Common configuration

A `[silico]` section contains setup for both `simuscop` and `varben` for clinvar variants insertion :

```toml
[silico]
# VCF containing clinvar variants. If not set, the VCF will be downloaded from NCBI
clinvar = "data/exp_raw/clinvar_col6a1.vcf.gz"
# Number of random clinvar variants to insert into the BAM file. Default is 1000
nb_variants = 2
# URL to fasta, or link to local version, otherwise download it from NCBI
# fasta =  "https://github.com/nf-core/test-datasets/blob/sarek3/data/genomics/homo_sapiens/genome/chr21/sequence/genome.fasta"
```

and a BAM file for varben is also needd in the common section as it can be used by simuscop or vaben, but for different purposes.

```toml
# Capture kit name (the BED must be defined in [capture] below)
capture = "agilent"
# A BAM file is always required for varben as it will will modify it
# Simuscop will define a sequencing profile from it if it does not exist.
# Can be a URL or a local path.
bam_file = "data/exp_raw/HG002.hiseq4000.wes-agilent.50x.dedup.grch38_nohardclip.bam"
```


### Simuscop specific configuration

```toml
# This section enables simuscop FASTQ generation (remove section to disable)
[silico.simuscop]
# Pre-built seqToProfile profile directory. If set, seqToProfile is skipped.
profile = "data/exp_raw/hiseq400-agilent-50x.profile"
# VCF of germline variants called from bam_file (e.g. via GATK HaplotypeCaller).
# Required when profile is absent; seqToProfile is run to build the profile.
#vcf = "data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
# Sequencing coverage
coverage = 50
```

Simuscop will generate a FASTQ according to a profile. `sherloxome` ships several pre-built profiles.

| Profile path                                | Sequencer  | Kit     | Depth |
| --                                          | --         | --      | --    |
| data/exp_raw/hiseq4000-agilent-50x.profile  | Hiseq 4000 | Agilent | 50x   |
| data/exp_raw/hiseq4000-idt-50x.profile      | Hiseq 4000 | IDT     | 50x   |
| data/exp_raw/hiseq4000-truseq-50x.profile   | Hiseq 4000 | Truseq  | 50x   |
| data/exp_raw/hiseq4000-agilent-75x.profile  | Hiseq 4000 | Agilent | 75x   |
| data/exp_raw/hiseq4000-idt-75x.profile      | Hiseq 4000 | IDT     | 75x   |
| data/exp_raw/hiseq4000-truseq-75x.profile   | Hiseq 4000 | Truseq  | 75x   |
| data/exp_raw/hiseq4000-agilent-100x.profile | Hiseq 4000 | Agilent | 100x  |
| data/exp_raw/hiseq4000-idt-100x.profile     | Hiseq 4000 | IDT     | 100x  |
| data/exp_raw/hiseq4000-truseq-100x.profile  | Hiseq 4000 | Truseq  | 100x  |
| data/exp_raw/novaseq-agilent-50x.profile    | novaseq    | Agilent | 50x   |
| data/exp_raw/novaseq-idt-50x.profile        | novaseq    | IDT     | 50x   |
| data/exp_raw/novaseq-truseq-50x.profile     | novaseq    | Truseq  | 50x   |
| data/exp_raw/novaseq-agilent-75x.profile    | novaseq    | Agilent | 75x   |
| data/exp_raw/novaseq-idt-75x.profile        | novaseq    | IDT     | 75x   |
| data/exp_raw/novaseq-truseq-75x.profile     | novaseq    | Truseq  | 75x   |
| data/exp_raw/novaseq-agilent-100x.profile   | novaseq    | Agilent | 100x  |
| data/exp_raw/novaseq-idt-100x.profile       | novaseq    | IDT     | 100x  |
| data/exp_raw/novaseq-truseq-100x.profile    | novaseq    | Truseq  | 100x  |


To create your own profile, a BAM, VCF are required

```toml
[silico.simuscop]
# VCF of germline variants called from bam_file (e.g. via GATK HaplotypeCaller).
# Required when profile is absent; seqToProfile is run to build the profile.
vcf = "data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
```

FASTQ will be generated in `data/exp_raw/simuscop_$CONFIG` as `$CONFIG_1.fq` and `$CONFIG_2.fq`

In the example above, the relevant part of the samplesheet is
```csv
silico-simuscop,agilent-col6a1_simuscop,1,data/exp_raw/simuscop_agilent-col6a1/agilent-col6a1_1.fq.gz,data/exp_raw/simuscop_agilent-col6a1/agilent-col6a1_2.fq.gz
```

Simuscop also requires a file with SNP. We download dbSNP and keep relevant variants [as documented here](021-dbsnp.md).

### Varben specific configuration

```toml
# This section enables varben BAM editing (remove section to disable)
[silico.varben]
mindepth = 30
```

FASTQ will be generated in `data/exp_raw/$PATIENT_$SEQUENCER_CONFIG`
In the example above, the relevant part of the samplesheet is
```csv
silico-varben,HG002_hiseq4000_agilent-col6a1_50x_nohardclip_varben,1,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_nohardclip_1.fq.gz,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_nohardclip_2.fq.gz
```

#### Varben algorithm

1. The BAM is downloaded (if `bam` is an URL) or verified locally
2. Hard-clipped reads are removed with `samtools` + `awk`
3. `nb_variants` ClinVar pathogenic SNVs are sampled from within the capture BED, enforcing ≥50 bp spacing between selected variants
4. `muteditor` inserts the variants into the BAM with random allele fractions (0.4–0.6)
5. The edited BAM is converted to paired FASTQ via `samtools fastq`
6. A VCF of successfully inserted variants is written alongside the FASTQ

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
