#!/bin/bash
#SBATCH --job-name=smk_herit_conda
#SBATCH --output=snakemake_herit_conda_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

SNAKEMAKE=/scratch.global/GDC/envs/snakemake/bin/snakemake
CONDA_PREFIX=/scratch.global/coffm049/tmp/conda_envs

$SNAKEMAKE \
    --directory /scratch.global/GDC/GDCGenomicsQC \
    --snakefile /scratch.global/GDC/GDCGenomicsQC/workflow/Snakefile \
    --configfile /scratch.global/GDC/GDCGenomicsQC/config/sandboxABCD.yaml \
    --use-conda \
    --conda-frontend mamba \
    --conda-prefix $CONDA_PREFIX \
    --cores 8 \
    --jobs 1 \
    --latency-wait 60 \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --until run_snpHerit \
    --config \
        conda-frontend=mamba \
        local-storage-prefix=/scratch.global/coffm049/tmp \
    2>&1 | tee snakemake_herit_conda_full.log
