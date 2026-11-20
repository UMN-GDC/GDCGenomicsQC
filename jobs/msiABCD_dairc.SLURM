#!/bin/bash
#SBATCH --job-name=dairc_lai_test
#SBATCH --output=dairc_lai_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use "/scratch.global/coffm049/GDCGenomicsQC/envs"
module load gdcgenomicsMSI

# Full pipeline through local ancestry (phasing + RFMIX, test mode):
gdcgenomicsqc --until run_lai --configfile "/scratch.global/coffm049/GDCGenomicsQC/config/msiABCD_dairc.yaml"

# TRACTOR run (local-ancestry-aware GWAS; requires LAI from above):
#   1. set gwas.pheno and add "tractor" to gwas.methods in the config
#   2. gdcgenomicsqc --until run_snpHerit --configfile "/scratch.global/coffm049/GDCGenomicsQC/config/msiABCD_dairc.yaml"
