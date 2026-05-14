#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module use /scratch.global/GDC/GDCGenomicsQC/envs
module load gdcgenomicsSandbox

gdcgenomicsqc /scratch.global/GDC/r25outputs/ABCD/full/f1.ldpruned.unrelated.grm.bin /scratch.global/GDC/r25outputs/ABCD/full/internal_pca_plink2.eigenvec --configfile ../config/sandboxABCD.yaml
