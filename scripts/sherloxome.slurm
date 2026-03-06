#!/bin/bash -l
# Fichier submission.SBATCH

#SBATCH --job-name="sherloxome-idt"
#SBATCH --output=%x.%J.out   ## %x=job name, %J=job id
#SBATCH --error=%x.%J.out
 # walltime (hh:mm::ss) max is 8 days
#SBATCH -t 24:00:00
#SBATCH --partition=smp
#SBATCH -c 1  ## request 16 cores (MAX is 32)
#SBATCH --mem=12G ## (MAX is 96G)
#SBATCH --mail-user=apraga@chu-besancon.fr
#SBATCH --mail-type=END,FAIL   # notify when job end/fail

# module load nix/2.11.0

# Otherwise job fails as it cannot write to $HOME/.nextflow
export NXF_HOME=/Work/Users/apraga/.nextflow

# user.name must be forced (again... our fix does not seem to work in nix)
nextflow run nf-core/sarek --input samplesheet.csv  -r 3.5.1 --outdir giab -c conf/slurm.config \
	--tools haplotypecaller --skip_tools haplotypecaller_filter	 \
	--wes --intervals ./capture/xgen-exome-hyb-panel-v2-targets-hg38.bed -bg \
	-resume	grave_spence
