#!/bin/bash -l        
#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=10g
#SBATCH --tmp=10g
#SBATCH -p agsmall
## SBATCH --mail-type=ALL  
## SBATCH --mail-user=x500@umn.edu 

module load plink
module load python3/3.9.3_anaconda2021.11_mamba
module load R/4.0.4
# module load R/4.2.2-openblas

# path_to_data=/home/gdc/shared/GDC_pipeline/K_branch/data/toySim/
# filepreffix=toy
path_to_data=/home/gdc/shared/SMILES/GSA/
filepreffix=SMILES_GSA
path_to_store_outputs=/home/gdc/shared/GDC_pipeline/

# Prefix for 
./src/main.sh --FILE ${filepreffix} --OLDBLD N --UPDTSO N --PATHTODATA ${path_to_data} --PATHTOSTOREOUTPUTS ${path_to_store_outputs} #Assumes they don't need their reference build updated and don't need alleles flipped. 

#./src/ancestryPrediction.sh #Output is ${filepreffix}_11 ? <-- doesn't work appropriately


./src/make_tables.sh --FILE ${filepreffix} --PATHTOSTOREOUTPUTS ${path_to_store_outputs} 

### Would be a good idea to also move the files in necessary_files.txt into the same data folder
## ./src/move_necessary_files.sh -- ............



# Generates PDF QC_report 
mkdir -p ${path_to_store_outputs}/results
Rscript ./src/QC_report.R ${path_to_store_outputs} ${filepreffix} 
# Errors here when switching working locations. Might be better to just point to the data rather than moving the working directory... 


## Currently not working so continuing with QC_report.R for testing purposes
# ./src/generate_report.sh # calls quarto render QCReport.qmd #For rendering the report

# Cleaning up repository
./src/repo_cleaner.sh --FILE ./data/default_list_of_files_to_remove.txt