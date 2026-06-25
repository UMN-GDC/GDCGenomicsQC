#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=10:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

# Snakemake auto-resolves the DAG from this single target —
# it builds QC → PCA → classification → ancestry QC → simulation → heritability:
gdcgenomicsqc run_snpHeritSimulated --configfile ../config/sandboxToy.yaml
