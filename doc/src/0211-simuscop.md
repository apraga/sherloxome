# Silico FASTQ with simuscop

```toml
# This section enables simuscop FASTQ generation (remove section to disable)
[silico.simuscop]
# Pre-built seqToProfile profile directory. The profile must follow the filenaming scheme
profile = "data/ref/hiseq4000_agilent_50x.profile"
# Sequencing coverage
coverage = 50
```

**Warning** : simuscop use a maximum coverage. The configuration above will be converted to an estimation for a mean coverage by dividing by 0.65 (empirical value).

Simuscop generates a FASTQ from a pre-built seqToProfile `profile`. `sherloxome` ships:

| Profile path                                | Sequencer  | Kit     | Depth |
| ------------------------------------------- | ---------- | ------- | ----- |
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

Variants inserted in the FASTQ are available as a VCF in `data/exp_raw`, for example
`data/exp_raw/nopatient_hiseq4000_agilent_50x_simuscop.vcf.gz`.

Simuscop will also add background SNPs from dbSNP on top of the ClinVar variants above.
See [the relevant section](022-dbsnp.md) for more information.

## GIAB example

`data/ref/hiseq4000_agilent_50x.profile` was itself built with `seqToProfile` from a real GIAB
HG002 BAM (HiSeq 4000, Agilent, 50x) :


```toml
[real]
patients   = ["HG002"]
sequencers = ["hiseq4000"]
captures   = ["agilent"]
depths     = [50]

[silico]
capture = "agilent"

[silico.simuscop]
profile = "data/ref/hiseq4000_agilent_50x.profile"
coverage = 50
```

`sherloxome setup` writes both rows into `samplesheet-agilent.csv`:

```csv
patient,sample,lane,fastq_1,fastq_2
HG002,HG002_hiseq4000_agilent_50x,1,https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes_agilent.50x.R1.fastq.gz,https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes_agilent.50x.R2.fastq.gz
silico-simuscop,nopatient_hiseq4000_agilent_50x_simuscop,1,data/exp_raw/nopatient_hiseq4000_agilent_50x_simuscop_1.fq.gz,data/exp_raw/nopatient_hiseq4000_agilent_50x_simuscop_2.fq.gz
```

Both rows can then go through the same sarek run and be benchmarked against the GIAB truth VCF
with `sherloxome benchmark` — see [Evaluate performance](05-benchmark.md).

## Building a new profile

`sherloxome` no longer builds seqToProfile profiles itself — `profile` must already exist on
disk. To add support for a new sequencer/capture/depth combination, run `seqToProfile` directly
(bundled in the `simuscop` Nix package, see `pkgs/simuscop`) against a real BAM and its called
variants.

For example, here's how it was done for GIAB HG002 on an Hiseq4000 with a depth of 50x and agilent capture kit :

```bash
base="https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/"
wget -nc "${base}/bam/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam.bai"  "${base}/bam/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam  "${base}/vcf/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes-agilent.50x.gatk4.grch38.vcf.gz"
gunzip  HG002.hiseq4000.wes-agilent.50x.gatk4.grch38.vcf.gz
```
Then run seqToProfile

```bash
seqToProfile  HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam -v HG002.hiseq4000.wes-agilent.50x.gatk4.grch38.vcf -t data/capture/agilent.targets.grch38.chr21.bed data/ref/hiseq4000_agilent_50x.profile

```

