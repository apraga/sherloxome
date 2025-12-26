#!/bin/bash

error() {
  echo "$@" 1>&2
}

fail() {
  error "$@"
  exit 1
}

# Compute the capture filename from the VCF full path
capture_bed() {
    if [[ $1 =~ agilent ]]; then
        echo "./capture/Agilent_SureSelect_All_Exons_v8_hg38_Regions.bed"
    else
        fail "No capture found in VCF filename $1"
    fi
}

# Compute the patient from the VCF full path
patient_from_vcf() {
    vcf=$(basename $1)
    for n in $(seq 1 7); do
        if [[ $vcf =~ "HG00$n" ]]; then
            echo "HG00$n"
            exit 0
        fi
    done
    fail "No patient found in VCF filename $1"
}

# Compute the truth filename from the patient name
truth_vcf() {
    matches=(./ref/*"$1"*.vcf.gz)
    nb=${#matches[@]}
    if [[ $nb == 1 ]]; then
        echo "${matches[0]}"
    else
        fail "Found $nb reference VCF for patient $1"
    fi
}

# Compute the truth filename from the patient name
truth_bed() {
    matches=(./ref/*"$1"*.bed)
    nb=${#matches[@]}
    if [[ $nb == 1 ]]; then
        echo "${matches[0]}"
    else
        fail "Found $nb reference bed for patient $1"
    fi
}


# Argument : input_VCF fasta output_dir
# Other parameters are deduced from the filename
happy() {
    patient=$(patient_from_vcf $1) || exit
    echo "Patient: $patient"
    truth_vcf=$(truth_vcf $patient) || exit
    truth_bed=$(truth_bed $patient) || exit
    capture_bed=$(capture_bed $1) || exit
	hap.py $truth_vcf $1 \
    -r $2 \
	--false-positives $truth_bed \
	--target-regions $capture_bed \
	--engine=vcfeval --engine-vcfeval-template ./sdf --threads 1 -o $3
}

# happy lol.vcf hg002-out
fasta=/Work/Groups/bisonex/dgenomes/genome-human/GCA_000001405.15_GRCh38_full_analysis_set.fna
vcf=./hg002/variant_calling/haplotypecaller/HG002-hiseq4000-agilent-50x/HG002-hiseq4000-agilent-50x.haplotypecaller.vcf.gz
# vcf=./HG002.hiseq4000.wes-agilent.50x.gatk4.grch38.vcf.gz
happy $vcf $fasta hg002-out
