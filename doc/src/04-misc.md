# Misc scripts

## Compute the coverage of multiple bam

`scripts/coverage` runs `mosdepth` (needs to be installed) on all .bam files in a directory and output the stats in an output directory.
Usage : 

    cargo run -p coverage --  -d BAMDIRECTORY -o OUTDIRECTORY

Ex: 

    cargo run -p coverage - -- -d ../../../bisonex/out/preprocessing/markduplicates/ -o bisonex-giab

Example of a slurm file (needs `cargo build --release` to be run beforehand)

```slurm
#!/bin/bash -l
# Fichier submission.SBATCH

#SBATCH --job-name="mosdepth"
#SBATCH --output=%x.%J.out   ## %x=job name, %J=job id
#SBATCH --error=%x.%J.out
 # walltime (hh:mm::ss) max is 8 days
#SBATCH -t 4:00:00
#SBATCH --partition=smp
#SBATCH -c 6  ## request 16 cores (MAX is 32)
#SBATCH --mem=12G ## (MAX is 96G)

cargo run -p coverage - -d ../../../bisonex/out/preprocessing/markduplicates/ -o bisonex-giab
```
