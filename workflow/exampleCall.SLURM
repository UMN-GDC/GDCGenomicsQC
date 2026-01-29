#!/bin/bash
#SBATCH --job-name=SM_Master
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time=07:00:00  # This is your 420 minutes

# Now call Snakemake inside this allocated space
snakemake --executor slurm --jobs 100 --profile my_slurm_profile
