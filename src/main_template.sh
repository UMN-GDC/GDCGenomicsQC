#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH -p msismall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=x500@umn.edu 
#SBATCH -o FLE.out
#SBATCH -e FLE.err
#SBATCH --job-name FLE

# This pipeline assumes the input is in plink binary format.

#################################### Specifying paths #########################################

# Hard-code the path to the Reference folder (containing reference dataset, other bash scripts, and programs' executables like CrossMap, GenomeHarmonizer, PRIMUS, and fraposa)
REF=/home/gdc/shared/GDC_pipeline/Ref
path_to_repo=PRPO
FILE=PND
NAME=FLE
WORK=WK
crossmap=CRSMP
genome_harmonizer=GNHRM
rfmix_option=RFMX
report_writer=RPT
custom_qc=CSTQC

cd ${WORK}

#################################################################################################

source /home/faird/shared/code/external/envs/miniconda3/load_miniconda3.sh
conda activate GDC_pipeline
module load plink
module load perl

############## Updating genome build and conducting strand alignment/allele flipping #############
#### Skipping everything until resume place when choosing to skip Crossmap ####
if [ ${crossmap} -eq 1 ]; then
  echo "(Step 1) Matching data to NIH's GRCh38 genome build"
  ${path_to_repo}/src/run_crossmap.sh ${WORK} ${REF} ${FILE} ${NAME} ${path_to_repo}
  file_to_use=study.${NAME}.lifted
  #plink --file ${file_to_use} --make-bed --out ${file_to_use} #Unsure if this is necessary
else  # Default behavior
  file_to_use=${FILE}/${NAME} #Original file
fi

#### Actual resume place for skipping updating genome build ####
# Break the dataset by chromosomes for faster processing in the next step (genome harmonizer)
if [ ${genome_harmonizer} -eq 1 ]; then
  echo "Begin genome harmonization"
  ${path_to_repo}/src/run_genome_harmonizer.sh ${WORK} ${REF} ${NAME} ${path_to_repo} ${file_to_use} #file_to_use is the primary change
  file_to_submit=$WORK/aligned/study.$NAME.lifted.aligned
else # Default behavior
  if [ ${crossmap} -eq 1 ]; then
    file_to_submit=study.$NAME.lifted #For using crossmap but not genome harmonizer
  else # Not using crossmap or genome harmonizer
    file_to_submit=${FILE}/${NAME} #Original file
  fi
fi
#######################################################################################################


###################################### QC #############################################################
echo "(Step 2) Standard variants and samples filtering"
# Run standard_QC.job with the appropriate parameters (full path to dataset name + output folder name)
cd $WORK
DATATYPE=full
if [ ${custom_qc} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  sbatch --wait ${WORK}/custom_qc.SLURM ${file_to_submit} ${DATATYPE} ${path_to_repo}
else # Default behavior
  sbatch --wait ${path_to_repo}/src/standard_QC.job ${file_to_submit} ${DATATYPE} ${path_to_repo}
fi
########################################################################################################


