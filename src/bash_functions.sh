#!/bin/bash


run_crossmap_if_needed() {
  local crossmap_check="$1"
  local path_to_repo="$2"
  local WORK="$3"
  local REF="$4"
  local FILE="$5"
  local NAME="$6"

  if [ ! -f "${crossmap_check}" ]; then
    echo "(Step 1) Matching data to NIH's GRCh38 genome build"
    sbatch --wait ${path_to_repo}/src/run_crossmap.sh ${WORK} ${REF} ${FILE} ${NAME} ${path_to_repo}
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
    sbatch --wait ${path_to_repo}/src/run_genome_harmonizer.sh ${WORK} ${REF} ${NAME} ${path_to_repo} ${file_to_use} 
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
    sbatch --wait ${path_to_repo}/src/initial_QC.sh ${file_to_submit} ${path_to_repo}
  fi
}

# run_standard_qc_if_needed ${file_to_check_qc} ${path_to_repo} ${file_to_submit} ${DATATYPE} ## Sample call


initial_qc_check_after_call() {
  local file_to_check_qc="$1"

  if [ ! -f "${file_to_check_qc}" ]; then
    echo "Standard QC steps have failed please check the error logs."
    exit 1
  fi
}

# standard_qc_check_after_call() ${file_to_check_qc} ## Sample call

run_standard_qc_if_needed() {
  local file_to_check_qc="$1"
  local WORK="$2"
  local NAME="$3"
  local path_to_repo="$4"
  local DATATYPE="$5"

  if [ ! -f "${file_to_check_qc}" ]; then
    sbatch --wait ${path_to_repo}/src/standard_QC.job ${WORK} ${NAME} ${path_to_repo} ${DATATYPE}
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
  local DATATYPE="$6"

  if [ ! -f "${king_check}" ]; then
    echo "(Step 3) Relatedness check"
    ${path_to_repo}/src/run_king.sh ${WORK} ${REF} ${NAME} ${path_to_repo} ${DATATYPE}
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
    sbatch --wait "${path_to_repo}/src/run_phase.sh" "${WORK}" "${REF}" "${NAME}" "${DATATYPE}" "${path_to_repo}"
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
    sbatch --wait "${path_to_repo}/src/run_rfmix.sh" "${WORK}" "${REF}" "${NAME}" "${path_to_repo}"
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
    sbatch --wait ${path_to_repo}/src/run_subpops.sh ${WORK} ${REF} ${NAME} ${path_to_repo}
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

  for DATATYPE in ${ETHNICS}; do
    mkdir -p $DATATYPE
    plink --bfile ${WORK}/aligned/study.${NAME}.lifted.aligned --keep ${WORK}/PCA/${DATATYPE} --make-bed --out ${WORK}/${DATATYPE}/study.${NAME}.${DATATYPE}.lifted.aligned
    if [ ${custom_qc} -eq 1 ]; then
    ## Will follow a pre-determined naming such as ${WORK}/custom_qc.SLURM
      sbatch ${WORK}/custom_qc.SLURM ${WORK}/aligned/study.${NAME}.${DATATYPE}.lifted.aligned ${DATATYPE} ${path_to_repo}
    else # Default behavior      
      sbatch ${path_to_repo}/src/per_ancestry_QC.job ${WORK} ${WORK}/aligned/study.${NAME}.${DATATYPE}.lifted.aligned ${path_to_repo}
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
  cp ${WORK}/aligned/*harmonizer*.txt ${WORK}/full/
  primus_file=$(find ${WORK} -type f -name "full.QC8.genome_network1.genome")
  cp ${primus_file} ${WORK}/full/primus_file.genome

  #3. move other directories into a temporary location called 'temp'
  # aligned, lifted, logs, PCA, relatedness, relatedness_OLD
  mkdir ${WORK}/temp
  mv -f ${WORK}/aligned ${WORK}/temp/
  mv -f ${WORK}/lifted ${WORK}/temp/
  mv -f ${WORK}/logs ${WORK}/temp/
  mv -f ${WORK}/Initial_QC ${WORK}/temp/Initial_QC
  mv -f ${WORK}/phased ${WORK}/temp/
  mv -f ${WORK}/rfmix ${WORK}/temp/
  mv -f ${WORK}/PCA ${WORK}/temp/
  mv -f ${WORK}/GAP_plots ${WORK}/temp/
  mv -f ${WORK}/LAP_plots ${WORK}/temp/
  mv -f ${WORK}/relatedness ${WORK}/temp/
  mv -f ${WORK}/relatedness_OLD ${WORK}/temp/
  mv -f ${WORK}/ancestry_estimation ${WORK}/temp/
  mv -f ${WORK}/*.out ${WORK}/temp/logs/out/
  mv -f ${WORK}/*.err ${WORK}/temp/logs/errors/

  mv ${WORK}/*.lifted* ${WORK}/temp/lifted #To clean up the working directory of unnecessary files 
}

# restructure_and_clean_outputs ${WORK} ${NAME} ## Sample call


