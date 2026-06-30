#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsLasso

#gdcgenomicsqc --until run_classifyAncestry --configfile "$HOME/GDCGenomicsQC/config/sandboxABCD_smokescreen.yaml"
#gdcgenomicsqc --until run_rfmix --configfile "$HOME/GDCGenomicsQC/config/sandboxABCD_smokescreen.yaml"
#gdcgenomicsqc --until run_ancestryQC --configfile "$HOME/GDCGenomicsQC/config/sandboxABCD_smokescreen.yaml"
gdcgenomicsqc --until run_snpHerit --configfile "$HOME/GDCGenomicsQC/config/lassoABCD_smokescreen.yaml"
