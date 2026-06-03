#!/bin/bash
#SBATCH --job-name=gdc_qc
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=7-00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

cd $SLURM_SUBMIT_DIR/GDCGenomicsQC/workflow

CONFIG="${CONFIG:-../config/config.yaml}"

snakemake --profile=../profiles/hpc \
    --configfile "$CONFIG" \
    --jobs 20
