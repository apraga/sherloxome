# Possibles errors

## Hap.py error about vcfeval
If hap.py raises an error with vcfeval about awk/hostname not found and --threads issues, force --threads to be less than 10. In SLURM, the number of threads can be set to the number of cores (so more than the maximum allowed by vcfeval).
This should not happy as the code forces vcfeval to run sequentially (for now).
