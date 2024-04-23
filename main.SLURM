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


# Start QC steps
./src/main.sh --FILE ${filepreffix} --OLDBLD N --UPDTSO N --PATHTODATA ${path_to_data} --PATHTOSTOREOUTPUTS ${path_to_store_outputs} #Assumes they don't need their reference build updated and don't need alleles flipped. 

#./src/ancestryPrediction.sh #Output is ${filepreffix}_11 ? <-- doesn't work appropriately

# Makes txt files to be used in the report based on log files
./src/make_tables.sh --FILE ${filepreffix} --PATHTOSTOREOUTPUTS ${path_to_store_outputs} 

# Moving specific log files that are used in the report into the proper final location
./src/move_necessary_files.sh --FILE ./data/default_list_of_files_to_move.txt --END_LOCATION ${path_to_store_outputs}/logs/



mkdir -p ${path_to_store_outputs}/results
#Rscript ./src/QC_report.R ${path_to_store_outputs} ${filepreffix}  
module load R #Need a newer version of R than 4.0.4 for quarto

## Makes automated pdf report
./src/generate_report.sh --FILE ${filepreffix} --PATHTOSTOREOUTPUTS ${path_to_store_outputs} # calls quarto render QCReport.qmd #For rendering the report

./src/move_necessary_files.sh --FILE ./data/default_list_of_files_to_move.txt --END_LOCATION ${path_to_store_outputs}/logs/

# Cleaning up repository /working directory
./src/repo_cleaner.sh --FILE ./data/default_list_of_files_to_remove.txt
