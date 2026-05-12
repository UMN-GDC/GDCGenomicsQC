#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=10:00:00

module use /scratch.global/GDC/GDCGenomicsQC/envs
module load gdcgenomicsSandbox

#gdcgenomicsqc run_classifyAncestry --configfile ../config/sandboxToy_noPredict.yaml
# gdcgenomicsqc run_ancestryQC --configfile ../config/sandboxToy_noPredict.yaml
#gdcgenomicsqc run_simulatePhenotype --configfile ../config/sandboxToy_noPredict.yaml
gdcgenomicsqc runAllEnabledPRS --configfile ../config/sandboxToy_noPredict.yaml --rerun-incomplete
