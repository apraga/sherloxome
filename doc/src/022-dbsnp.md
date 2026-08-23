# dbSNP setup for simuscop

To add background SNPs to generated FASTQ by simuscop, we use a trimmed-down version of
dbSNP, keeping only
1. variants in the selected capture kit
2. frequent variants with AF > 1% (dbSNP's own `COMMON` flag: minor allele frequency >= 1%
   in at least one 1000Genomes population, with 2+ founders contributing)
3. absent from our clinvar selection, so the two variant lists never collide on the same
   position

The CLI requires an existing file for filtered dbSNP data to match the capture kit and only keep common variants.
dbSNP data for the agilent, IDT and truseq capture kit are already alvailable in the source code. For another capture kit, run

```bash
  # 0. Download dbSNP
  curl "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
 # 1. capture bed: chr -> RefSeq (for `bcftools view -R`)
 awk -F'\t' 'NR==FNR{gsub(/\r/,""); m[$1]=$2; next} $1 in m{print m[$1]"\t"$2"\t"$3}' \
     "$MAPPING" "$BED" > ${OUT}.regions.bed

 # 2. build the reverse map: RefSeq -> chr (for `bcftools annotate --rename-chrs`)
 awk -F'\t' '{gsub(/\r/,""); print $2"\t"$1}' "$MAPPING" > ${OUT}.rename.txt

 # 3. filter dbSNP to the capture kit, common variants, SNVs only
 bcftools view -v snps -i 'COMMON=1' -R ${OUT}.regions.bed "$DBSNP" -Oz -o ${OUT}.refseq.vcf.gz

 # 4. convert result back to chr notation
 bcftools annotate --rename-chrs ${OUT}.rename.txt ${OUT}.refseq.vcf.gz -Oz -o ${OUT}.vcf.gz

 # 5. index
 bcftools index -t ${OUT}.vcf.gz
```

Rename the output file to `data/exp_raw/dbsnp_{CAPTURE}_common.vcf.gz`.

The CLI will remove clinvar variants from this VCF and write the final result as a plain text file for Simuscop.
