#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00  # Give the master enough time to finish the whole pipe

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
#snakemake --profile=../profiles/hpc --configfile=../config/toy.yaml --rerun-incomplete run_pca
snakemake --profile=../profiles/hpc --configfile=../config/toy.yaml --rerun-incomplete run_classifyAncestry

