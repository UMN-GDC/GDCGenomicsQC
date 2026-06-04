#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=01:30:00
#SBATCH --job-name=buildsnakemake
#SBATCH --output=snakemake.out
#SBATCH --partition=compute-light

source ~/miniconda3/etc/profile.d/conda.sh

conda create -n snakemake_env -c conda-forge -c bioconda snakemake snakemake-executor-plugin-slurm
