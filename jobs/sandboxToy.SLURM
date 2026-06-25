#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=10:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

# Full pipeline (uncomment when ready):
# gdcgenomicsqc run_classifyAncestry --configfile ../config/sandboxToy.yaml
# gdcgenomicsqc run_ancestryQC run_relatedness run_PCAinternal run_simulatePhenotype run_snpHeritSimulated --configfile ../config/sandboxToy.yaml

# Testing MASH on simulated data (comment out after testing):
gdcgenomicsqc run_snpHeritSimulated --configfile ../config/sandboxToy.yaml
