#!/bin/bash
#SBATCH --job-name=mash_estimate
#SBATCH --output=mash_%j.log
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=2880

module use "$HOME/GDCGenomicsQC/envs"
module load mashMSI

mash estimate --argfile {ARGFILE}