#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=4G
#SBATCH --time=72:00:00

# Edit this path to match your conda installation
source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

# The magic flag is --executor slurm
snakemake --profile=../profiles/hpc
