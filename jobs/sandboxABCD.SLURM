#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

gdcgenomicsqc --until run_classifyAncestry --configfile ../config/sandboxABCD.yaml
#gdcgenomicsqc --until run_snpHerit --configfile ../config/sandboxABCD2.yaml
