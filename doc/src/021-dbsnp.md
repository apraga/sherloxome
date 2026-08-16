# dbSNP setup for simuscop

To add SNPs to generated fastq by simuscop, we use a trimmed-down version of dbSNP, keeping only
1. variants in the selected capture kit
2. frequent variants with AF > 1%
3. absent from our clinvar selection

Note that  dbSNP use Refseq naming convention for chromosome, so it requires a mapping file. This file is available in data/ref/chromosome_mapping_GRCh38.p41.txt. It can be generated with

```bash
curl -s "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt" | awk -F'\t' '$2=="assembled-molecule"{print $NF"\t"$7}' > data/ref/chromosome_mapping_GRCh38.p14.txt
```
