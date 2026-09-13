#!/bin/bash
set -euo pipefail

base="https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38"

paths=(
    "bam/hiseq4000/wes_agilent/50x/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam.bai"
    "bam/hiseq4000/wes_idt/100x/HG002.hiseq4000.wes_idt.100x.dedup.bam.bai"
    "bam/hiseq4000/wes_idt/50x/HG002.hiseq4000.wes_idt.50x.dedup.bam.bai"
    "bam/hiseq4000/wes_idt/75x/HG002.hiseq4000.wes_idt.75x.dedup.bam.bai"
    "bam/hiseq4000/wes_truseq/50x/HG002.hiseq4000.wes-truseq.50x.dedup.grch38.bam.bai"
    "bam/hiseq4000/wes_truseq/75x/HG002.hiseq4000.wes-truseq.75x.dedup.grch38.bam.bai"
    "bam/novaseq/wes_agilent/100x/HG002.novaseq.wes-agilent.100x.dedup.grch38.bam.bai"
    "bam/novaseq/wes_agilent/50x/HG002.novaseq.wes-agilent.50x.dedup.grch38.bam.bai"
    "bam/novaseq/wes_agilent/75x/HG002.novaseq.wes-agilent.75x.dedup.grch38.bam.bai"
    "bam/novaseq/wes_idt/100x/HG002.novaseq.wes_idt.100x.dedup.bam.bai"
    "bam/novaseq/wes_idt/50x/HG002.novaseq.wes_idt.50x.dedup.bam.bai"
    "bam/novaseq/wes_idt/75x/HG002.novaseq.wes_idt.75x.dedup.bam.bai"
    "bam/novaseq/wes_truseq/100x/HG002.novaseq.wes_truseq.100x.dedup.bam.bai"
    "bam/novaseq/wes_truseq/50x/HG002.novaseq.wes_truseq.50x.dedup.bam.bai"
    "bam/novaseq/wes_truseq/75x/HG002.novaseq.wes_truseq.75x.dedup.bam.bai"
)

for path in "${paths[@]}"; do
    wget -nc "${base}/${path}"
done
