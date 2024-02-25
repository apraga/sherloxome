# TODO Evaluate performance

## For reference patients
In `scripts/compare` will get all the HG***.vcf.gz in a directory, match them against reference patients using gold vcf and bed dataset.
It also needs capture kits in `baid2020/bed`

Usage : 

    cargo run -p compare -- -d VCFDIRECTORY -o OUTDIRECTORY

Ex: 

    cargo run -p compare -- -d ../../../bisonex/out/preprocessing/markduplicates/ -o bisonex-giab

Example of a slurm file (needs `cargo build --release` to be run beforehand)

```slurm

#!/bin/bash -l
# Fichier submission.SBATCH

#SBATCH --job-name="compare-exome-baid2020"
#SBATCH --output=%x.%J.out   ## %x=job name, %J=job id
#SBATCH --error=%x.%J.out
 # walltime (hh:mm::ss) max is 8 days
#SBATCH -t 4:00:00
#SBATCH --partition=smp
#SBATCH -c 6  ## request 16 cores (MAX is 32)
#SBATCH --mem=12G ## (MAX is 96G)

module load nix/2.11.0
cargo run -p compare --  -d ../../baid2020/grch38/vcf -o giab-baid2020
```

## TODO For synthetic patient
## TODO For synthetic data

