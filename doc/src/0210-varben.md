# FASTQ from patient with Varben

 This section enables varben BAM editing (remove section to disable)
```toml
[silico.varben]
bam_files = ["data/exp_raw/HG002_hiseq4000_agilent_50x.bam", "data/exp_raw/HG002_novaseq_agilent_75x.bam"]
mindepth = 30
```

`bam_files` takes one entry per patient/sequencer/depth combination you want to edit — all
must be for the same capture kit as `[silico] capture`, since that's what determines the
sampled clinvar panel shared across them (see [Common configuration](021-silico.md)). To
cover several capture kits, run `setup` once per kit with a different `[silico] capture` /
`bam_files` pair.

FASTQ will be generated in `data/exp_raw/$PATIENT_$SEQUENCER_CONFIG`, one pair per BAM.
In the example above, the relevant part of the samplesheet is
```csv
silico-varben,HG002_hiseq4000_agilent-col6a1_50x_nohardclip_varben,1,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_nohardclip_1.fq.gz,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_nohardclip_2.fq.gz
```

**Some variants may not in the FASTQ**. Those which failed to be inserted are available in `data/exp_raw` as txt file. For example  `data/exp_raw/HG002_hiseq4000_agilent_50x_varben_failed.txt`.
Variants successfully inserted in the BAM are available as a VCF is `data/exp_raw`, for example `data/exp_raw/HG002_hiseq4000_agilent_50x_varben.vcf.gz`. 

## Varben algorithm

The clinvar panel (`nb_variants` ClinVar pathogenic SNVs sampled from within the capture BED,
enforcing ≥50 bp spacing) is sampled once and reused identically for every BAM in
`bam_files` — this keeps runs comparable across patients/sequencers/depths, since any
difference in recovery is then attributable to the run's conditions rather than to a
different sample of variants. For each BAM in `bam_files`:

1. The BAM's filename is parsed to recover its patient/sequencer/depth (see below)
2. Hard-clipped reads are removed with `samtools` + `awk`
3. `muteditor` inserts the (shared) variant panel into the BAM with random allele fractions (0.4–0.6)
4. The edited BAM is converted to paired FASTQ via `samtools fastq`
5. A VCF of successfully inserted variants is written alongside the FASTQ

Each BAM in `bam_files` must already exist locally and be renamed to follow our
[filenaming scheme](050-filenaming.md) — URLs are not supported, since the filename itself
is how sherloxome recovers the run's patient/sequencer/depth. BAM files can be found
[on Google Cloud for data from Baid et al, 2020](https://console.cloud.google.com/storage/browser/brain-genomics-public/research/sequencing/grch38/bam;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false). For example, [HG002 data sequenced on Hiseq 4000 with Agilent capture kit](https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam) should be downloaded and renamed with :

```bash
mv HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam HG002_hiseq4000_agilent_50x.bam
```

`scripts/download-bam.rs` and `scripts/rename_to_underscore.sh` automate downloading and
renaming several such BAMs at once.

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


