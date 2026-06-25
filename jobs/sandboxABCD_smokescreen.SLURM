#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

#gdcgenomicsqc --until run_classifyAncestry --configfile ../config/sandboxABCD_smokescreen.yaml
#gdcgenomicsqc --until run_rfmix --configfile ../config/sandboxABCD_smokescreen.yaml
gdcgenomicsqc --until run_classifyAncestry run_ancestryQC run_snpHerit --configfile ../config/sandboxABCD_smokescreen.yaml
