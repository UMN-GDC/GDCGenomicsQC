#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=48:00:00
#SBATCH -p msismall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=x500@umn.edu 
#SBATCH -o FLE.out
#SBATCH -e FLE.err
#SBATCH --job-name FLE

# This pipeline assumes the input is in plink binary format.

#################################### Specifying paths #################################################

# Hard-code the path to the Reference folder (containing reference dataset, other bash scripts, and programs' executables like CrossMap, GenomeHarmonizer, PRIMUS, and fraposa)
REF=/home/gdc/public/Ref
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

####################################### Environment Setup ##############################################

conda activate /home/gdc/public/envs/GDC_pipeline
source ${path_to_repo}/src/bash_functions.sh # Helper functions
module load plink
module load perl

############## Updating genome build and conducting strand alignment/allele flipping ###################
if [ ${crossmap} -eq 1 ]; then
  file_to_use=study.${NAME}.lifted
  crossmap_check=${WORK}/${file_to_use}.bim
  run_crossmap_if_needed ${crossmap_check} ${path_to_repo} ${WORK} ${REF} ${FILE} ${NAME}
  crossmap_check_after_call ${crossmap_check}

  else  # Default behavior
  file_to_use=${FILE}/${NAME} #Original file
fi
############## Genome harmonizer section 
if [ ${genome_harmonizer} -eq 1 ]; then
  file_to_submit=$WORK/aligned/study.$NAME.lifted.aligned
  run_genome_harmonizer_if_needed ${file_to_submit} ${path_to_repo} ${WORK} ${REF} ${NAME} ${file_to_use}
  genome_harmonizer_check_after_call ${file_to_submit}
else # Default behavior
  if [ ${crossmap} -eq 1 ]; then
    file_to_submit=study.$NAME.lifted #For using crossmap but not genome harmonizer
  else # Not using crossmap or genome harmonizer
    file_to_submit=${FILE}/${NAME} #Original file
  fi
fi
#######################################################################################################


###################################### QC #############################################################
echo "Variants and samples filtering"
# Run standard_QC.job with the appropriate parameters (full path to dataset name + output folder name)
cd $WORK
DATATYPE=full
if [ ${custom_qc} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  sbatch --wait ${WORK}/custom_qc.SLURM ${file_to_submit} ${DATATYPE} ${path_to_repo}
else # Default behavior
  file_to_check_qc=${WORK}/${DATATYPE}/${DATATYPE}.QC8.bim
  run_standard_qc_if_needed ${file_to_check_qc} ${path_to_repo} ${file_to_submit} ${DATATYPE}
  standard_qc_check_after_call ${file_to_check_qc}
fi
########################################################################################################


######################################## Pedigree ######################################################
primus_check=$WORK/relatedness/study.$NAME.unrelated.bim
run_primus_if_needed ${primus_check} ${path_to_repo} ${WORK} ${REF} ${NAME} ${DATATYPE}
primus_check_after_call ${primus_check}
#########################################################################################################


######################################## Phasing ########################################################
echo "(Step 4) Phasing"
if [ ${rfmix_option} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  run_phasing_if_needed ${WORK} ${REF} ${NAME} ${path_to_repo} 
  phasing_check_after_call ${WORK} ${NAME}
else
  echo "Skip phasing and move to Fraposa"
fi
#########################################################################################################


######################################## Ethnicity ######################################################
echo "(Step 5) ancestry estimate"
if [ ${rfmix_option} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  run_rfmix_if_needed ${WORK} ${REF} ${NAME} ${path_to_repo}
  rfmix_check_after_call ${WORK}
else # Alternative behavior
  ${path_to_repo}/src/run_fraposa.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
fi
##########################################################################################################


###################################### Subpopulations ####################################################
echo "(Step 6) Subpopulations"
subpop_check=${WORK}/PCA/study.${NAME}.unrelated.comm.popu
if [ ${rfmix_option} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  run_subpopulations_if_needed ${subpop_check} ${path_to_repo} ${WORK} ${REF} ${NAME}
else # Alternative behavior
  echo "Skip subpopulations"
fi
subpop_check_after_call ${subpop_check}
##########################################################################################################


############################################ PCA #########################################################
echo "PCA"
sbatch --wait ${path_to_repo}/src/run_pca.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
##########################################################################################################


######################################## Ancestry Plots ##################################################
echo "(Step 7) ancestry plots"
if [ ${rfmix_option} -eq 1 ]; then
  sbatch --wait ${path_to_repo}/src/run_rfmix_plots.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
  rm -r ${WORK}/visualization
  Rscript ${path_to_repo}/src/plot_pca.R ${WORK}
else # Alternative behavior
  echo "Plot module only for rfmix"
fi
#########################################################################################################


################### Subset data based on Ethnicity and Rerun QC (Step 2) on the subsets #################
cd ${WORK}
if [ ${rfmix_option} -eq 1 ]; then
  ETHNICS=$(awk '{print $3}' ${WORK}/PCA/study.${NAME}.unrelated.comm.popu | sort | uniq)
else # Alternative behavior
  ETHNICS=$(awk -F'\t' '{print $3}' ${WORK}/PCA/study.${NAME}.unrelated.comm.popu | sort | uniq)
fi

subset_ancestries_run_standard_qc "${ETHNICS}" ${WORK} ${NAME} ${custom_qc} ${path_to_repo}
##Putting in to wait until the jobs are done
wait_for_ancestry_qc_to_finish
###########################################################################################################

  
########################## Restructuring and cleaning up for the report writer ############################
restructure_and_clean_outputs ${WORK} ${NAME}

#4. execute run_generate_reports.sh ##
module load R/4.4.0-openblas-rocky8
if [ ${report_writer} -eq 1 ]; then
  ${path_to_repo}/src/run_generate_reports.sh ${WORK} ${path_to_repo}
fi

