# dbSNP setup for simuscop

To add background SNPs to generated FASTQ by simuscop, we use a trimmed-down version of
dbSNP, keeping only
1. variants in the selected capture kit
2. frequent variants with AF > 1% (dbSNP's own `COMMON` flag: minor allele frequency >= 1%
   in at least one 1000Genomes population, with 2+ founders contributing)
3. absent from our clinvar selection, so the two variant lists never collide on the same
   position

Rather than downloading the ~30GB dbSNP VCF, `sherloxome` runs a remote `tabix` query
restricted to the capture kit regions, which only fetches the required byte ranges over HTTP.

Note that dbSNP uses Refseq naming convention for chromosome, so a mapping file is required
to translate to/from `chr` notation. This file is available in
`data/ref/chromosome_mapping_GRCh38.p14.txt`. It can be (re)generated with

```bash
curl -s "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt" | awk -F'\t' '$2=="assembled-molecule"{print $NF"\t"$7}' > data/ref/chromosome_mapping_GRCh38.p14.txt
```

## Configuration

Enable this by adding a `[silico.simuscop.dbsnp]` section:

```toml
[silico.simuscop.dbsnp]
# Local bgzipped+tabix-indexed dbSNP VCF, or a remote URL.
# Defaults to the latest NCBI GRCh38 dbSNP release.
# vcf = "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
# chr -> Refseq accession mapping. Defaults to data/ref/chromosome_mapping_GRCh38.p14.txt
# mapping = "data/ref/chromosome_mapping_GRCh38.p14.txt"
```

## Output files

- `data/exp_raw/dbsnp_{capture}.vcf.gz` (+ `.tbi`) — the filtered dbSNP variants, in `chr`
  notation, sorted and bgzipped
- `data/exp_raw/dbsnp_{capture}.snp` — the same variants converted to simuscop's own SNP
  format and passed as the `snp` parameter to `simuReads` (see `SimuSCoP_User_Guide.pdf` §4.2.2)

Both are skipped (existing files reused) on subsequent runs, matching the ClinVar sampling
step.
