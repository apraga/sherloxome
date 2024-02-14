
# Run a pipeline for reference patients

## TODO Generate samplesheet

## TODO run pipeline

## Compare output

In `scripts/compare-exome` will get all the HG***.vcf.gz in a directory, match them against reference patients using gold vcf and bed dataset.
It also needs capture kits in `baid2020/bed`

Usage : 

    compare-exome -d VCFDIRECTORY -o OUTDIRECTORY

Ex: 

    cargo run -- -d ../../../bisonex/out/preprocessing/markduplicates/ -o bisonex-giab

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
./target/release/compare-exome -d ../../baid2020/grch38/vcf -o giab-baid2020
```

 # Misc scripts

## Compute the coverage of multiple bam

scripts/coverage subproject run mosdepth (needs to be installed) on all .bam files in a directory and output the stats in an output directory.
Usage : 

    coverage -d BAMDIRECTORY -o OUTDIRECTORY

Ex: 

    cargo run -- -d ../../../bisonex/out/preprocessing/markduplicates/ -o bisonex-giab

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

./target/release/coverage -d ../../../bisonex/out/preprocessing/markduplicates/ -o bisonex-giab
```

# Possibles errors

If hap.py raises an error with vcfeval about awk/hostname not found and --threads issues, force --threads to be less than 10. In SLURM, the number of threads can be set to the number of cores (so more than the maximum allowed by vcfeval).
This should not happy as the code forces vcfeval to run sequentially (for now).
