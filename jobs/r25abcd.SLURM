#!/bin/bash
#SBATCH --job-name=smk_master
#SBATCH --output=snakemake_%j.log
#SBATCH --mem=1G
#SBATCH --time=72:00:00

module purge
MODULEPATH=/scratch.global/coffm049/GDCGenomicsQC/envs:$MODULEPATH
module load gdcgenomicsqc

gdcgenomicsqc run_classifyAncestry --configfile ../config/r25abcd.yaml 
