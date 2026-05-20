#!/bin/bash
#SBATCH --job-name=smk_prs
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

gdcgenomicsqc -n --until runAllEnabledPRS --configfile ../config/minimal_prs_external.yaml