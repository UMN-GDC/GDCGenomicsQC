#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use /scratch.global/coffm049/GDCGenomicsQC/envs
module load gdcgenomicsSandbox

gdcgenomicsqc assembleRef --configfile ../config/sandboxABCD.yaml --force
