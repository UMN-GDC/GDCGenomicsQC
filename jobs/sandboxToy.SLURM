#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=10:00:00

module use "$HOME/GDCGenomicsQC/envs"
module load gdcgenomicsSandbox

# Single target — Snakemake auto-resolves the full DAG:
# QC → PCA → UMAP → classification → ancestry QC → simulation → heritability
#gdcgenomicsqc run_classifyAncestry run_snpHeritSimulated --configfile "$HOME/GDCGenomicsQC/config/sandboxToy.yaml"
gdcgenomicsqc run_classifyAncestry run_initialQC --configfile "$HOME/GDCGenomicsQC/config/sandboxToy2.yaml"
