# FASTQ from patient with Varben

 This section enables varben BAM editing (remove section to disable)
```toml
[silico.varben]
mindepth = 30
```

FASTQ will be generated in `data/exp_raw/$PATIENT_$SEQUENCER_CONFIG`
In the example above, the relevant part of the samplesheet is
```csv
silico-varben,HG002_hiseq4000_agilent-col6a1_50x_nohardclip_varben,1,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_nohardclip_1.fq.gz,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_nohardclip_2.fq.gz
```

**Some variants may not in the FASTQ**. Those which failed to be inserted are available in `data/exp_raw` as txt file. For example  `data/exp_raw/HG002_hiseq4000_agilent_50x_varben_failed.txt`.
Variants successfully inserted in the BAM are available as a VCF is `data/exp_raw`, for example `data/exp_raw/HG002_hiseq4000_agilent_50x_varben.vcf.gz`. 

## Varben algorithm

1. The BAM is downloaded (if `bam` is an URL) or verified locally
2. Hard-clipped reads are removed with `samtools` + `awk`
3. `nb_variants` ClinVar pathogenic SNVs are sampled from within the capture BED, enforcing ≥50 bp spacing between selected variants
4. `muteditor` inserts the variants into the BAM with random allele fractions (0.4–0.6)
5. The edited BAM is converted to paired FASTQ via `samtools fastq`
6. A VCF of successfully inserted variants is written alongside the FASTQ

BAM files can be found [on Google Cloud for data from Baid et al, 2020](https://console.cloud.google.com/storage/browser/brain-genomics-public/research/sequencing/grch38/bam;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false). Those must be renamed to follow our [filenaming scheme](050-filenaming.md). For example,  [HG002 data sequenced on Hiseq 4000 with Agilent capture kit](https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam) should be renamed with :

```bash
mv HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam HG002_hiseq4000_agilent_50x.bam
```

# Output files

By default, output is written to `data/exp_raw/`. Override with `outdir`:

```toml
[silico]
outdir = "silico"
```

| File                                    | Description                           |
|-----------------------------------------|---------------------------------------|
| `clinvar_{capture}.vcf.gz`              | Sampled ClinVar variants (truth VCF)  |
| `clinvar_{capture}.mut`                 | Variants to insert                    |
| `clinvar_{capture}_varben_failed.txt`          | Variants to insert                    |
| `varben/edit.sorted.bam`                | BAM with variants injected            |
| `{sample}.vcf.gz`                       | VCF of successfully inserted variants |
| `{sample}_1.fq.gz` / `{sample}_2.fq.gz` | Paired FASTQ ready for the pipeline   |


