#!/bin/bash

if [[ -n "$APPTAINER_NAME" || -n "$SINGULARITY_NAME" ]]; then
    module load apptainer
    SUBMIT_SUF="apptainer exec $APPTAINER_NAME"
else
    SUBMIT_SUF=""
fi


run_crossmap_if_needed() {
  local crossmap_check="$1"
  local path_to_repo="$2"
  local WORK="$3"
  local REF="$4"
  local FILE="$5"
  local NAME="$6"

  if [ ! -f "${crossmap_check}" ]; then
    echo "(Step 1) Matching data to NIH's GRCh38 genome build"
    sbatch --wait $SUBMIT_SUF ${path_to_repo}/src/run_crossmap.sh ${WORK} ${REF} ${FILE} ${NAME} ${path_to_repo}
  fi
}

# run_crossmap_if_needed "/path/to/checkfile" "/path/to/repo" "$WORK" "$REF" "$FILE" "$NAME" ## Sample call


crossmap_check_after_call() {
  local crossmap_check="$1"

  if [ ! -f "${crossmap_check}" ]; then
    echo "Crossmap has failed please check the error logs."
    exit 1
  fi
}

# crossmap_check_after_call "/path/to/checkfile" ## Sample call


run_genome_harmonizer_if_needed() {
  local file_to_submit="$1"
  local path_to_repo="$2"
  local WORK="$3"
  local REF="$4"
  local NAME="$5"
  local file_to_use="$6"

  if [ ! -f "${file_to_submit}.bim" ]; then
    echo "Begin genome harmonization"
    sbatch --wait $SUBMIT_SUF ${path_to_repo}/src/run_genome_harmonizer.sh ${WORK} ${REF} ${NAME} ${path_to_repo} ${file_to_use} 
  fi
}

# run_genome_harmonizer_if_needed ${file_to_submit} ${path_to_repo} "$WORK" "$REF" "$NAME" "$file_to_use" ## Sample call


genome_harmonizer_check_after_call() {
  local file_to_submit="$1"

  if [ ! -f "${file_to_submit}.bim" ]; then
    echo "Genome Harmonizer has failed please check the error logs."
    exit 1
  fi
}

# genome_harmonizer_check_after_call ${file_to_submit} ## Sample call 

run_initial_qc_if_needed() {
  local file_to_check_qc="$1"
  local path_to_repo="$2"
  local file_to_submit="$3"

  if [ ! -f "${file_to_check_qc}" ]; then
    sbatch --wait $SUBMIT_SUF ${path_to_repo}/src/initial_QC.sh ${file_to_submit} ${path_to_repo}
  fi
}

# run_standard_qc_if_needed ${file_to_check_qc} ${path_to_repo} ${file_to_submit} ${DATATYPE} ## Sample call


initial_qc_check_after_call() {
  local file_to_check_qc="$1"

  if [ ! -f "${file_to_check_qc}" ]; then
    echo "Initial QC steps have failed please check the error logs."
    exit 1
  fi
}

# standard_qc_check_after_call() ${file_to_check_qc} ## Sample call

run_standard_qc_if_needed() {
  local file_to_check_qc="$1"
  local WORK="$2"
  local NAME="$3"
  local REF="$4"
  local path_to_repo="$5"
  local DATATYPE="$6"
  local CHECK_SEX="$7"

  if [ ! -f "${file_to_check_qc}" ]; then
    sbatch --wait $SUBMIT_SUF ${path_to_repo}/src/standard_QC.job ${WORK} ${NAME} ${REF} ${path_to_repo} ${DATATYPE} ${CHECK_SEX}
  fi
}

# run_standard_qc_if_needed ${file_to_check_qc} ${path_to_repo} ${file_to_submit} ${DATATYPE} ## Sample call


standard_qc_check_after_call() {
  local file_to_check_qc="$1"

  if [ ! -f "${file_to_check_qc}" ]; then
    echo "Standard QC steps have failed please check the error logs."
    exit 1
  fi
}

# standard_qc_check_after_call() ${file_to_check_qc} ## Sample call


run_primus_if_needed() {
  local primus_check="$1"
  local path_to_repo="$2"
  local WORK="$3"
  local REF="$4"
  local NAME="$5"
  local DATATYPE="$6"

  if [ ! -f "${primus_check}" ]; then
    echo "(Step 3) Relatedness check"
    ${path_to_repo}/src/run_primus.sh ${WORK} ${REF} ${NAME} ${path_to_repo} ${DATATYPE}
  fi
}

# run_primus_if_needed ${primus_check} ${path_to_repo} ${WORK} ${REF} ${NAME} ${DATATYPE} ## Sample call


primus_check_after_call() {
  local primus_check="$1"

  if [ ! -f "${primus_check}" ]; then
    echo "Primus relatedness check has failed please check the error logs."
    exit 1
  fi
}

# primus_check_after_call ${primus_check} ## Sample call

run_king_if_needed() {
  local king_check="$1"
  local path_to_repo="$2"
  local WORK="$3"
  local REF="$4"
  local NAME="$5"
  local COMB="$6"

  if [ ! -f "${king_check}" ]; then
    echo "(Step 3) Relatedness check"
    ${path_to_repo}/src/run_king.sh ${WORK} ${REF} ${NAME} ${path_to_repo} ${COMB}
  fi
}

# run_king_if_needed ${king_check} ${path_to_repo} ${WORK} ${REF} ${NAME} ${DATATYPE} ## Sample call


king_check_after_call() {
  local king_check="$1"

  if [ ! -f "${king_check}" ]; then
    echo "King relatedness check has failed please check the error logs."
    exit 1
  fi
}

# king_check_after_call ${king_check} ## Sample call

run_pca_ir_if_needed() {
  local pcair_check="$1"
  local path_to_repo="$2"
  local WORK="$3"
  local REF="$4"
  local NAME="$5"

  if [ ! -f "${pcair_check}" ]; then
    echo "(Step 3) Running PC-AiR and PC-Relate"
    ${path_to_repo}/src/run_pcair.sh ${WORK} ${REF} ${NAME} ${path_to_repo} 
  fi
}

# run_pca_ir_if_needed ${pcair_check} ${path_to_repo} ${WORK} ${REF} ${NAME} ## Sample call

pca_ir_check_after_call() {
  local pcair_check="$1"

  if [ ! -f "${pcair_check}" ]; then
    echo "PC-AiR and PC-Relate relatedness check has failed please check the error logs."
    exit 1
  fi
}

# pca_ir_check_after_call ${pcair_check} ## Sample call


run_phasing_if_needed() {
  local WORK="$1"
  local REF="$2"
  local NAME="$3"
  local path_to_repo="$4"
  local DATATYPE="$5"

  # Generate the list of expected phase files
  local phase_files=()
  for CHR in {1..22}; do
    phase_files+=("${WORK}/phased/${NAME}.chr${CHR}.phased.vcf.gz")
  done

  # Check if all phase files exist
  local all_exist=true
  for phase_check in "${phase_files[@]}"; do
    if [ ! -f "$phase_check" ]; then
      all_exist=false
      break
    fi
  done

  # If any file is missing, run the phasing script
  if ! $all_exist; then
    sbatch --wait $SUBMIT_SUF "${path_to_repo}/src/run_phase.sh" "${WORK}" "${REF}" "${NAME}" "${DATATYPE}" "${path_to_repo}"
  fi
}

# run_phasing_if_needed ${WORK} ${REF} ${NAME} ${path_to_repo} ## Sample call


phasing_check_after_call() {
  local WORK="$1"
  local NAME="$2"

  # Generate the list of expected phase files
  local phase_files=()
  for CHR in {1..22}; do
    phase_files+=("${WORK}/phased/${NAME}.chr${CHR}.phased.vcf.gz")
  done

  for phase_check in "${phase_files[@]}"; do
    if [ ! -f "$phase_check" ]; then
      echo "Phasing has failed, please check the error logs."
      exit 1
    fi
  done
}

# phasing_check_after_call ${WORK} ${NAME} ## Sample call


run_rfmix_if_needed() {
  local WORK="$1"
  local REF="$2"
  local NAME="$3"
  local path_to_repo="$4"

  # Generate the list of expected RFMix files
  local rfmix_files=()
  for CHR in {1..22}; do
    rfmix_files+=("${WORK}/rfmix/ancestry_chr${CHR}.rfmix.Q")
  done

  # Check if all RFMix files exist
  local all_exist=true
  for rfmix_check in "${rfmix_files[@]}"; do
    if [ ! -f "$rfmix_check" ]; then
      all_exist=false
      break
    fi
  done

  # If any file is missing, run the RFMix script
  if ! $all_exist; then
    sbatch --wait $SUBMIT_SUF "${path_to_repo}/src/run_rfmix.sh" "${WORK}" "${REF}" "${NAME}" "${path_to_repo}"
  fi
}

# run_rfmix_if_needed ${WORK} ${REF} ${NAME} ${path_to_repo} ## Sample call


rfmix_check_after_call() {
  local WORK="$1"

  # Generate the list of expected RFMix files
  local rfmix_files=()
  for CHR in {1..22}; do
    rfmix_files+=("${WORK}/rfmix/ancestry_chr${CHR}.rfmix.Q")
  done

  for rfmix_check in "${rfmix_files[@]}"; do
    if [ ! -f "$rfmix_check" ]; then
      echo "Rfmix has failed, please check the error logs."
      exit 1
    fi
  done
}

# rfmix_check_after_call ${WORK} ## Sample call


run_subpopulations_if_needed() {
  local subpop_check="$1"
  local path_to_repo="$2"
  local WORK="$3"
  local REF="$4"
  local NAME="$5"

  if [ ! -f "${subpop_check}" ]; then
    sbatch --wait $SUBMIT_SUF ${path_to_repo}/src/run_subpops.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
  fi
}

# run_subpopulations_if_needed ${subpop_check} ${path_to_repo} ${WORK} ${REF} ${NAME} ## Sample call


subpop_check_after_call() {
  local subpop_check="$1"

  if [ ! -f "${subpop_check}" ]; then
    echo "Subpopulations estimation has failed please check the error logs."
    exit 1
  fi
}

# subpop_check_after_call ${subpop_check} ## Sample call


subset_ancestries_run_standard_qc() {
  local ETHNICS="$1"
  local WORK="$2"
  local NAME="$3"
  local custom_qc="$4"
  local path_to_repo="$5"
  local ORIG=${6:-NULL}
  for DATATYPE in ${ETHNICS}; do
    mkdir -p $DATATYPE

    if [[ -f "${WORK}/aligned/study.${NAME}.lifted.aligned.bed" ]]; then
      plink --bfile ${WORK}/aligned/study.${NAME}.lifted.aligned --keep ${WORK}/PCA/${DATATYPE} --make-bed --out ${WORK}/${DATATYPE}/study.${NAME}.lifted.aligned
      if [ ${custom_qc} -eq 1 ]; then
      ## Will follow a pre-determined naming such as ${WORK}/custom_qc.SLURM
        sbatch ${WORK}/custom_qc.SLURM ${WORK}/${DATATYPE}/study.${NAME}.${DATATYPE}.lifted.aligned ${DATATYPE} ${path_to_repo}
      else # Default behavior      
        sbatch $SUBMIT_SUF ${path_to_repo}/src/per_ancestry_QC.job ${WORK}/${DATATYPE}/study.${NAME}.${DATATYPE}.lifted.aligned ${DATATYPE} ${path_to_repo}
      fi
    else 
      plink --bfile $ORIG --keep ${WORK}/PCA/${DATATYPE} --make-bed --out ${WORK}/${DATATYPE}/study.${NAME}
      if [ ${custom_qc} -eq 1 ]; then
      ## Will follow a pre-determined naming such as ${WORK}/custom_qc.SLURM
        sbatch ${WORK}/custom_qc.SLURM ${WORK}/${DATATYPE}/study.${NAME}.${DATATYPE} ${DATATYPE} ${path_to_repo}
      else # Default behavior      
        sbatch $SUBMIT_SUF ${path_to_repo}/src/per_ancestry_QC.job ${WORK}/${DATATYPE}/study.${NAME} ${DATATYPE} ${REF} ${path_to_repo}
      fi
    fi
  done
}

# subset_ancestries_run_standard_qc ${ETHNICS} ${WORK} ${NAME} ${custom_qc} ${path_to_repo} ## Sample call


wait_for_ancestry_qc_to_finish() {
  local jobs_remaining=$(squeue --me | grep -c QC)
  echo "${jobs_remaining} jobs remaining at the start of this waiting loop"
  
  local x=1
  while [ ${jobs_remaining} -gt 0 ]
  do
    sleep 1m
    jobs_remaining=$(squeue --me | grep -c QC)
    echo "${jobs_remaining} after waiting for ${x} minutes"
    ((x++))
  done
}

# wait_for_ancestry_qc_to_finish ## Sample call


restructure_and_clean_outputs() {
  local WORK="$1"
  local NAME="$2"

  #1. move over png and .popu file from PCA directory into the 'full' directory
  cp ${WORK}/PCA/study.${NAME}*popu ${WORK}/full/
  cp ${WORK}/PCA/*png ${WORK}/full/

  #2. move the genome_harmonizer_full_log.txt into the 'full' directory
  harmonizer_file=$(find ${WORK} -name "*harmonizer*.txt" | head -n 1)
  cp -u ${harmonizer_file} ${WORK}/full/
  king_file=$(find ${WORK} -type f -name "kinships.kin0" | head -n 1)
  cp -u ${king_file} ${WORK}/full/kinships.kin0

  #3. move other directories into a temporary location called 'temp'
  # aligned, lifted, logs, PCA, relatedness, relatedness_OLD
  mkdir ${WORK}/temp
  mv -f ${WORK}/Initial_QC ${WORK}/temp/Initial_QC
  if [[ -d "${WORK}/aligned" ]]; then
    mv -f ${WORK}/aligned ${WORK}/temp/
  fi 
  if [[ -d "${WORK}/lifted" ]]; then
    mv -f ${WORK}/lifted ${WORK}/temp/
  fi
  if [[ -d "${WORK}/logs" ]]; then
    mv -f ${WORK}/logs ${WORK}/temp/
  fi
  if [[ -d "${WORK}/phased" ]]; then
    mv -f ${WORK}/phased ${WORK}/temp/
  fi
  if [[ -d "${WORK}/rfmix" ]]; then
    mv -f ${WORK}/rfmix ${WORK}/temp/
  fi
  if [[ -d "${WORK}/PCA" ]]; then
    mv -f ${WORK}/PCA ${WORK}/temp/
  fi
  if [[ -d "${WORK}/gds" ]]; then
    mv -f ${WORK}/gds ${WORK}/temp/
  fi
  if [[ -d "${WORK}/visualization" ]]; then
    mv -f ${WORK}/visualization ${WORK}/temp/
  fi
  if [[ -d "${WORK}/pca_ir" ]]; then
    mv -f ${WORK}/pca_ir ${WORK}/temp/
  fi
  if [[ -d "${WORK}/GAP_plots" ]]; then
    mv -f ${WORK}/GAP_plots ${WORK}/temp/
  fi
  if [[ -d "${WORK}/LAP_plots" ]]; then
    mv -f ${WORK}/LAP_plots ${WORK}/temp/
  fi
  if [[ -d "${WORK}/relatedness" ]]; then
    mv -f ${WORK}/relatedness ${WORK}/temp/
  fi
  if [[ -d "${WORK}/relatedness_OLD" ]]; then
    mv -f ${WORK}/relatedness_OLD ${WORK}/temp/
  fi
  if [[ -d "${WORK}/ancestry_estimation" ]]; then
    mv -f ${WORK}/ancestry_estimation ${WORK}/temp/
  fi
  mkdir -p ${WORK}/temp/logs/out/
  mkdir -p ${WORK}/temp/logs/errors/
  mv -f ${WORK}/*.out ${WORK}/temp/logs/out/
  mv -f ${WORK}/*.err ${WORK}/temp/logs/errors/

  mv ${WORK}/*.lifted* ${WORK}/temp/lifted #To clean up the working directory of unnecessary files 
}

# restructure_and_clean_outputs ${WORK} ${NAME} ## Sample call


