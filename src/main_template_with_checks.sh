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
  
  crossmap_check=${WORK}/${file_to_use}.bim
  if [ ! -f "${crossmap_check}" ]; then
    echo "Crossmap has failed please check the error logs."
    exit 1
  fi
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
  
  if [ ! -f "${file_to_submit}.bim" ]; then
    echo "Genome Harmonizer has failed please check the error logs."
    exit 1
  fi
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
  
  file_to_check_qc=${WORK}/${DATATYPE}/${DATATYPE}.QC8.bim
  if [ ! -f "${file_to_check_qc}" ]; then
    echo "Standard QC steps have failed please check the error logs."
    exit 1
  fi
fi
########################################################################################################


######################################## Pedigree ######################################################
echo "(Step 3) Relatedness check"
${path_to_repo}/src/run_primus.sh ${WORK} ${REF} ${NAME} ${path_to_repo} ${DATATYPE}

primus_check=$WORK/relatedness/study.$NAME.unrelated.bim
if [ ! -f "${primus_check}" ]; then
  echo "Primus relatedness check has failed please check the error logs."
  exit 1
fi
#########################################################################################################


######################################## Ethnicity ######################################################
echo "(Step 4) PCA"
if [ ${rfmix_option} -eq 1 ]; then
  ## requires a text file that has all of the flags and specifications
  sbatch --wait ${path_to_repo}/src/run_rfmix.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
else # Alternative behavior
  ${path_to_repo}/src/run_fraposa.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
fi

pca_check=${WORK}/PCA/study.${NAME}.unrelated.comm.popu
if [ ! -f "${pca_check}" ]; then
  echo "PCA software has failed please check the error logs."
  exit 1
fi
#########################################################################################################


################### Subset data based on Ethnicity and Rerun QC (Step 2) on the subsets #################
cd ${WORK}
if [ ${rfmix_option} -eq 1 ]; then
  ETHNICS=$(awk '{print $3}' ${WORK}/PCA/study.${NAME}.unrelated.comm.popu | sort | uniq)
else # Alternative behavior
  ETHNICS=$(awk -F'\t' '{print $3}' ${WORK}/PCA/study.${NAME}.unrelated.comm.popu | sort | uniq)
fi

for DATATYPE in ${ETHNICS}; do
    plink --bfile ${WORK}/aligned/study.${NAME}.lifted.aligned --keep ${WORK}/PCA/${DATATYPE} --make-bed --out ${WORK}/aligned/study.${NAME}.${DATATYPE}.lifted.aligned
    if [ ${custom_qc} -eq 1 ]; then
    ## Will follow a pre-determined naming such as ${WORK}/custom_qc.SLURM
      sbatch ${WORK}/custom_qc.SLURM ${WORK}/aligned/study.${NAME}.${DATATYPE}.lifted.aligned ${DATATYPE} ${path_to_repo}
    else # Default behavior
      sbatch ${path_to_repo}/src/standard_QC.job ${WORK}/aligned/study.${NAME}.${DATATYPE}.lifted.aligned ${DATATYPE} ${path_to_repo}
    fi
done
###########################################################################################################

##Putting in to wait until the jobs are done
jobs_remaining=$(squeue --me | grep QC | wc -l)
echo "${jobs_remaining} jobs remaining at the start of this waiting loop"
x=1
while [ ${jobs_remaining} -gt 0 ]
do
  sleep 1m
  jobs_remaining=$(squeue --me | grep QC | wc -l)
  echo "${jobs_remaining} after waiting for ${x} minutes"
  ((x++))
done
  
########################## Restructuring and cleaning up for the report writer ############################
#1. move over png and .popu file from PCA directory into the 'full' directory
cp ${WORK}/PCA/study.${NAME}*popu ${WORK}/full/
cp ${WORK}/PCA/study.${NAME}*png ${WORK}/full/

#2. move the genome_harmonizer_full_log.txt into the 'full' directory
cp ${WORK}/aligned/*harmonizer*.txt ${WORK}/full/
primus_file=$(find ${WORK} -type f -name "full.QC8_cleaned.genome")
cp -v ${primus_file} ${WORK}/full/primus_file.genome

#3. move other directories into a temporary location called 'temp'
# aligned, lifted, logs, PCA, relatedness, relatedness_OLD
mkdir ${WORK}/temp
mv -f ${WORK}/aligned ${WORK}/temp/
mv -f ${WORK}/lifted ${WORK}/temp/
mv -f ${WORK}/logs ${WORK}/temp/
mv -f ${WORK}/phased ${WORK}/temp/
mv -f ${WORK}/rfmix ${WORK}/temp/
mv -f ${WORK}/PCA ${WORK}/temp/
mv -f ${WORK}/relatedness ${WORK}/temp/
mv -f ${WORK}/relatedness_OLD ${WORK}/temp/
mv -f ${WORK}/*.out ${WORK}/logs/out/
mv -f ${WORK}/*.err ${WORK}/logs/errors/

mv ${WORK}/*.lifted* ${WORK}/temp #To clean up the working directory of unnecessary files 

#4. execute run_generate_reports.sh ## Need to make this optional ##
module load R/4.4.0-openblas-rocky8
if [ ${report_writer} -eq 1 ]; then
  ${path_to_repo}/src/run_generate_reports.sh ${WORK} ${path_to_repo}
fi
