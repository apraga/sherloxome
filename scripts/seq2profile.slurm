#!/bin/bash -l
# Fichier submission.SBATCH

#SBATCH --job-name="seqToProfile"
#SBATCH --output=%x.%J.out   ## %x=job name, %J=job id
#SBATCH --error=%x.%J.out
 # walltime (hh:mm::ss) max is 8 days
#SBATCH -t 4:00:00
#SBATCH --partition=smp
#SBATCH -c 1  ## request 16 cores (MAX is 32)
#SBATCH --mem=16G ## (MAX is 96G)
#SBATCH --mail-user=apraga@chu-besancon.fr
#SBATCH --mail-type=END,FAIL   # notify when job end/fail

module load nix
./result/bin/seqToProfile  -b ../bisonex/out/preprocessing/recalibrated/2300346867_63118093-NA12878-GRCh38/2300346867_63118093-NA12878-GRCh38.recal.bam -t ../bisonex/capture/Twist_Exome_Core_Covered_Targets_hg38.bed -v ../bisonex/out/call_variant/haplotypecaller/2300346867_63118093-NA12878-GRCh38/2300346867_63118093-NA12878-GRCh38.haplotypecaller.vcf -r genome/GCA_000001405.15_GRCh38_full_analysis_set.fna > na12878-bisonex-grch38.model
