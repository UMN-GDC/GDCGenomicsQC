#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=4G
#SBATCH --time=72:00:00  # Give the master enough time to finish the whole pipe

module use /scratch.global/coffm049/GDCGenomicsQC/envs
module load gdcgenomicsMSI

# The magic flag is --executor slurm (or --cluster for older versions)
snakemake --profile=../profiles/hpc run_rfmix
