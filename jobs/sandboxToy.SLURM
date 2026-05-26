#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=10:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

gdcgenomicsqc run_classifyAncestry --configfile ../config/sandboxToy.yaml --container-image-prep
# gdcgenomicsqc run_ancestryQC --configfile ../config/sandboxToy.yaml 
#gdcgenomicsqc run_simulatePhenotype --configfile ../config/sandboxToy.yaml
# gdcgenomicsqc runAllEnabledPRS --configfile ../config/sandboxToy.yaml --rerun-incomplete 
