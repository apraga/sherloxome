# Silico FASTQ with simuscop]

``` toml
 This section enables simuscop FASTQ generation (remove section to disable)
[silico.simuscop]
 Pre-built seqToProfile profile directory. If set, seqToProfile is skipped.
profile = "data/exp_raw/hiseq400-agilent-50x.profile"
 VCF of germline variants called from bam_file (e.g. via GATK HaplotypeCaller).
 Required when profile is absent; seqToProfile is run to build the profile.
vcf = "data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
 Sequencing coverage 
coverage = 50
```

**Warning** : simuscop use a maximum coverage. The configuration above will be converted to an estimation for a mean coverage by dividing by 0.65 (empirical value).

Simuscop will generate a FASTQ according to a profile. `sherloxome` ships several pre-built profiles.

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
`data/exp_raw/nopatient_hiseq4000_agilent_50x_simuscop.vcf.gz`

## Creating your own simuscop profile 

To create your own profile, a BAM, VCF are required

``` toml
[silico.simuscop]
 VCF of germline variants called from bam_file (e.g. via GATK HaplotypeCaller).
 Required when profile is absent; seqToProfile is run to build the profile.
vcf = "data/ref/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
```

FASTQ will be generated in `data/exp_raw/simuscop_$CONFIG` as `$CONFIG_1.fq` and `$CONFIG_2.fq`

In the example above, the relevant part of the samplesheet is

``` csv
silico-simuscop,agilent-col6a1_simuscop,1,data/exp_raw/simuscop_agilent-col6a1/agilent-col6a1_1.fq.gz,data/exp_raw/simuscop_agilent-col6a1/agilent-col6a1_2.fq.gz
```

Simuscop will also add background SNPs from dbSNP on top of the ClinVar variants above. It requires dbSNP data filtered See [the relevant section](022-dbsnp) for more information.
