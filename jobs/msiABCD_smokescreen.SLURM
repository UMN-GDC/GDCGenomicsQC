#!/bin/bash
#SBATCH --job-name=smoke_master
#SBATCH --output=smoke_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use "/scratch.global/coffm049/GDCGenomicsQC/envs"
module load gdcgenomicsMSI

#gdcgenomicsqc --until run_classifyAncestry --configfile ../config/sandboxABCD_imputed.yaml
#gdcgenomicsqc --until run_rfmix --configfile ../config/sandboxABCD_imputed.yaml
#gdcgenomicsqc --until run_ancestryQC --configfile ../config/sandboxABCD_imputed.yaml
gdcgenomicsqc --until run_snpHerit --configfile ../config/msiABCD_smokescreen2.yaml
