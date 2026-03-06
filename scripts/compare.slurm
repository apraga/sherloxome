#!/bin/bash -l
# Fichier submission.SBATCH

#SBATCH --job-name="sherloxome-agilent-hiseq400"
#SBATCH --output=%x.%J.out   ## %x=job name, %J=job id
#SBATCH --error=%x.%J.out
 # walltime (hh:mm::ss) max is 8 days
#SBATCH -t 24:00:00
#SBATCH --partition=smp
#SBATCH -c 1  ## request 16 cores (MAX is 32)
#SBATCH --mem=12G ## (MAX is 96G)
#SBATCH --mail-user=apraga@chu-besancon.fr
#SBATCH --mail-type=END,FAIL   # notify when job end/fail

module load nix/2.11.0

make compare
