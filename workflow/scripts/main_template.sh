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
REF=/projects/standard/gdc/public/Ref
path_to_repo=PRPO
FILE=PND
NAME=FLE
WORK=WK
crossmap=CRSMP
genome_harmonizer=GNHRM
king=0
rfmix_option=RFMX
report_writer=RPT
custom_qc=CSTQC
combine_related=COMB
custom_ancestry=CSTANC
CHECK_SEX=CHECK_SEXX

cd ${WORK}

# Check environment if we are currently running inside Apptainer
if [[ -n "$APPTAINER_NAME" || -n "$SINGULARITY_NAME" ]]; then
    echo "Environment: Inside Container ($APPTAINER_NAME)"
    path_to_github_repo=/app/GDCGenomicsQC
    EXE="apptainer exec $APPTAINER_NAME"
else
    echo "Environment: Native Host"
    EXE="apptainer exec $APPTAINER_NAME"
fi




####################################### Environment Setup ##############################################
source /projects/standard/gdc/public/envs/load_miniconda3-2.sh
conda activate gdcPipeline
source ${path_to_repo}/src/bash_functions.sh # Helper functions

############## Updating genome build and conducting strand alignment/allele flipping ###################
if [ ${crossmap} -eq 1 ]; then
  file_to_use=study.${NAME}.lifted
  crossmap_check=${WORK}/${file_to_use}.bim
  run_crossmap_if_needed "${crossmap_check}" ${path_to_repo} ${WORK} ${REF} ${FILE} ${NAME}
  crossmap_check_after_call ${crossmap_check}

  else  # Default behavior
  file_to_use=${FILE}/${NAME} #Original file
fi
############## Genome harmonizer section 
if [ ${genome_harmonizer} -eq 1 ]; then
  file_to_submit=$WORK/aligned/study.$NAME.lifted.aligned
  run_genome_harmonizer_if_needed "${file_to_submit}" ${path_to_repo} ${WORK} ${REF} ${NAME} ${file_to_use}
  genome_harmonizer_check_after_call ${file_to_submit}
else # Default behavior
  if [ ${crossmap} -eq 1 ]; then
    file_to_submit=study.$NAME.lifted #For using crossmap but not genome harmonizer
  else # Not using crossmap or genome harmonizer
    file_to_submit=${FILE}/${NAME} #Original file
  fi
fi


#########################################################################################################

###################################### QC #############################################################

########################################################################################################


######################################## Phasing ########################################################
echo "(Step 5) Phasing"
if [ ${rfmix_option} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  run_phasing_if_needed ${WORK} ${REF} ${NAME} ${path_to_repo} ${DATATYPE}
  phasing_check_after_call ${WORK} ${NAME} ${DATATYPE}
else
  echo "Skip phasing and move to Fraposa"
fi
#########################################################################################################


######################################## Ethnicity ######################################################
echo "(Step 6) ancestry estimate"
if [ ${rfmix_option} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  run_rfmix_if_needed ${WORK} ${REF} ${NAME} ${path_to_repo}
  rfmix_check_after_call ${WORK}
else # Alternative behavior
  ${path_to_repo}/src/run_fraposa.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
fi
##########################################################################################################


######################################## Ancestry Plots ##################################################
echo "(Step 7) ancestry plots"
if [ ${rfmix_option} -eq 1 ]; then
  $EXE --wait ${path_to_repo}/src/run_rfmix_plots.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
  rm -r ${WORK}/visualization
else # Alternative behavior
  echo "Plot module only for rfmix"
fi
##########################################################################################################


############################################ PCA #########################################################
echo "(Step 8) PCA"
pcair_check=$WORK/pca_ir/${NAME}_pcaobj.RDS
run_pca_ir_if_needed ${pcair_check} ${path_to_repo} ${WORK} ${REF} ${NAME}
pca_ir_check_after_call ${pcair_check}
$EXE --wait ${path_to_repo}/src/run_pca.sh ${WORK} ${REF} ${NAME} ${path_to_repo} 
##########################################################################################################


###################################### Subpopulations ####################################################
echo "(Step 9) Subpopulations"
subpop_check=${WORK}/ancestry_estimation/study.${NAME}.unrelated.comm.popu
if [ ${rfmix_option} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  run_subpopulations_if_needed "${subpop_check}" ${path_to_repo} ${WORK} ${REF} ${NAME}
  Rscript ${path_to_repo}/src/plot_pca.R ${WORK}
else # Alternative behavior
  echo "Skip subpopulations"
fi
subpop_check_after_call ${subpop_check}
##########################################################################################################


# 
# 
# 
# ################### Subset data based on Ethnicity and Rerun QC (Step 2) on the subsets #################
# cd ${WORK}
# if [ ${rfmix_option} -eq 1 ]; then
#   ETHNICS=$(awk '{print $3}' ${WORK}/ancestry_estimation/study.${NAME}.unrelated.comm.popu | sort | uniq)
# else # Alternative behavior
#   ETHNICS=$(awk -F'\t' '{print $3}' ${WORK}/ancestry_estimation/study.${NAME}.unrelated.comm.popu | sort | uniq)
# fi
# 
# cp ${WORK}/ancestry_estimation/* ${WORK}/PCA/
# subset_ancestries_run_standard_qc "${ETHNICS}" ${WORK} ${NAME} ${custom_qc} ${path_to_repo} $FILE/$NAME
# ##Putting in to wait until the jobs are done
# wait_for_ancestry_qc_to_finish
# ###########################################################################################################

  
########################## Restructuring and cleaning up for the report writer ############################
restructure_and_clean_outputs ${WORK} ${NAME}
# 
# if [ ${report_writer} -eq 1 ]; then
#   ${path_to_repo}/src/run_generate_reports.sh ${WORK} ${path_to_repo}
# fi

