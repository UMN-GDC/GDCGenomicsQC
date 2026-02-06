#!/bin/bash
#SBATCH --job-name=SM_Master
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000       # doesn't take much memory for the scheduler script
#SBATCH --time=16:00:00  # This script needs to persist to schedule all downstream scripts

module load conda
conda activate snakemake
#snakemake --executor slurm --jobs 25 --profile my_slurm_profile

snakemake --cores=4 --use-conda \
    --configfile </path/to/config.yaml> --directory </path/to/GDCGenomicsQC/workflow> --snakefile </path/to/GDCGenomicsQC/workflow/Snakefile> \
    --executor slurm \
    --jobs 25
