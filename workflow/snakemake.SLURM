#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00  # Give the master enough time to finish the whole pipe

source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
snakemake --profile=../profiles/hpc --configfile=../config/abcd_config.yaml --rerun-incomplete run_classifyAncestry
