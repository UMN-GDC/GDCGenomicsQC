#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=10:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

gdcgenomicsqc run_ancestryQC run_snpHeritSimulated --configfile "$HOME/GDCGenomicsQC/config/sandboxToy.yaml"
