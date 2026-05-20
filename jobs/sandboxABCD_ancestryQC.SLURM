#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

gdcgenomicsqc --until run_ancestryQC --configfile ../config/sandboxABCD.yaml
