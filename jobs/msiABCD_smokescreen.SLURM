#!/bin/bash
#SBATCH --job-name=smoke_master
#SBATCH --output=smoke_%j.log
#SBATCH --mem=1G
#SBATCH --time=3:00:00

module use "/scratch.global/coffm049/GDCGenomicsQC/envs"
module load gdcgenomicsMSI

#gdcgenomicsqc --until run_classifyAncestry --configfile "/scratch.global/coffm049/GDCGenomicsQC/config/sandboxABCD_imputed.yaml"
#gdcgenomicsqc --until run_rfmix --configfile "/scratch.global/coffm049/GDCGenomicsQC/config/sandboxABCD_imputed.yaml"
#gdcgenomicsqc --until run_ancestryQC --configfile "/scratch.global/coffm049/GDCGenomicsQC/config/sandboxABCD_imputed.yaml"
gdcgenomicsqc --until run_snpHerit --configfile "/scratch.global/coffm049/GDCGenomicsQC/config/msiABCD_smokescreen.yaml"
