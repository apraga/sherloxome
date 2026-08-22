# FASTQ from patient with Varben

```toml
 This section enables varben BAM editing (remove section to disable)
[silico.varben]
mindepth = 30
```

FASTQ will be generated in `data/exp_raw/$PATIENT_$SEQUENCER_CONFIG`
In the example above, the relevant part of the samplesheet is
```csv
silico-varben,HG002_hiseq4000_agilent-col6a1_50x_nohardclip_varben,1,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_nohardclip_1.fq.gz,data/exp_raw/HG002_hiseq4000_agilent-col6a1_50x_nohardclip_2.fq.gz
```

## Varben algorithm

1. The BAM is downloaded (if `bam` is an URL) or verified locally
2. Hard-clipped reads are removed with `samtools` + `awk`
3. `nb_variants` ClinVar pathogenic SNVs are sampled from within the capture BED, enforcing ≥50 bp spacing between selected variants
4. `muteditor` inserts the variants into the BAM with random allele fractions (0.4–0.6)
5. The edited BAM is converted to paired FASTQ via `samtools fastq`
6. A VCF of successfully inserted variants is written alongside the FASTQ

# Output files

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

