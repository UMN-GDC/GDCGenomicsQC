#!/bin/bash -l        
#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH -p agsmall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=x500@umn.edu 

module load plink
module load python3/3.9.3_anaconda2021.11_mamba
module load R/4.2.2-openblas

# Prefix for 
./src/main.sh --FILE SMILES_GSA --OLDBLD N --UPDTSO N #Assumes they don't need their reference build updated and don't need alleles flipped. 
./src/ancestryPrediction.sh

./src/make_tables.sh --FILE SMILES_GSA 