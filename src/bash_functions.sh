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
    "${path_to_repo}/src/run_crossmap.sh" "${WORK}" "${REF}" "${FILE}" "${NAME}" "${path_to_repo}"
  fi
}

# run_crossmap_if_needed "/path/to/checkfile" "/path/to/repo" "$WORK" "$REF" "$FILE" "$NAME"

crossmap_check_after_call() {
  local crossmap_check="$1"

  if [ ! -f "${crossmap_check}" ]; then
    echo "Crossmap has failed please check the error logs."
    exit 1
  fi
}

# crossmap_check_after_call "/path/to/checkfile"

run_genome_harmonizer_if_needed() {
  local file_to_submit="$1"
  local path_to_repo="$2"
  local WORK="$3"
  local REF="$4"
  local FILE="$5"
  local NAME="$6"
  local file_to_use="$7"

  if [ ! -f "${file_to_submit}.bim" ]; then
    echo "Begin genome harmonization"
    ${path_to_repo}/src/run_genome_harmonizer.sh ${WORK} ${REF} ${NAME} ${path_to_repo} ${file_to_use} #file_to_use is the primary change
  fi
}

# run_genome_harmonizer_if_needed ${file_to_submit} ${path_to_repo} "$WORK" "$REF" "$FILE" "$NAME" "$file_to_use"
