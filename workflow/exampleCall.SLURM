#!/bin/bash
#SBATCH --job-name=SM_Master
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time=16:00:00  # This is your 420 minutes

# Now call Snakemake inside this allocated space
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
#snakemake --executor slurm --jobs 25 --profile my_slurm_profile
snakemake --executor slurm --jobs 25 --configfile=../config/config.yaml --use-conda run_train_gnomix
