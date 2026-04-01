#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00  # Give the master enough time to finish the whole pipe

source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

# The magic flag is --executor slurm (or --cluster for older versions)
snakemake --profile=../profiles/hpc --rerun-incomplete --set-resources "runtime=60"
