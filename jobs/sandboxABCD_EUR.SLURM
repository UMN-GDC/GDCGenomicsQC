#!/bin/bash
#SBATCH --job-name=smk_abcd_eur
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

# Runs QC → PCA → ancestry-specific QC for the EUR subset from self-report labels
gdcgenomicsqc --until run_ancestryQC --configfile ../config/sandboxABCD_EUR.yaml
