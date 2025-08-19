#!/bin/bash

# === Input Validation ===
if [ "$#" -ne 13 ]; then
  echo "Usage: $0 <input_dir> <input_file> <repo_path> <user_x500> <working_dir> <use_crossmap> <use_genome_harmonizer> <use_king> <use_rfmix> <make_report> <custom_qc> <custom_ancestry> <local_modules>"
  exit 1
fi

# === Assign Inputs ===
path_to_input_directory="$1"
input_file_name="$2"
path_github_repo="$3"
user_x500="$4"
desired_working_directory="$5"
using_crossmap="$6"
using_genome_harmonizer="$7"
using_king="$8"
using_rfmix="$9"
making_report="${10}"
custom_qc="${11}"
custom_ancestry="${12}"
local_modules="${13}"

# === Output Path Setup ===
output="${desired_working_directory}/${input_file_name}_wrapper.sh"

echo "[INFO] Creating output directory if it doesn't exist..."
mkdir -p "${desired_working_directory}" || {
  echo "[ERROR] Failed to create working directory: ${desired_working_directory}"
  exit 1
}

# === Template File Check ===
template="${path_github_repo}/src/main_template.sh"

if [ ! -f "${template}" ]; then
  echo "[ERROR] Template file not found: ${template}"
  exit 1
fi


output=${desired_working_directory}/${input_file_name}_wrapper.sh
mkdir -p ${desired_working_directory}


sed -e "s|PRPO|${path_github_repo}|g" \
    -e "s|PND|${path_to_input_directory}|g" \
    -e "s|FLE|${input_file_name}|g" \
    -e "s|x500|${user_x500}|g" \
    -e "s|WK|${desired_working_directory}|g" \
    -e "s|CRSMP|${using_crossmap}|g" \
    -e "s|GNHRM|${using_genome_harmonizer}|g" \
    -e "s|KING|${using_king}|g" \
    -e "s|RFMX|${using_rfmix}|g" \
    -e "s|RPT|${making_report}|g" \
    -e "s|CSTQC|${custom_qc}|g" \
    -e "s|CSTANC|${custom_ancestry}|g" \
    -e "s|LMDL|${local_modules}|g" \
    "${template}" > "${output}"

#cp -v ${path_github_repo}/src/main_template.sh ${output}
#sed -i 's@PND@'${path_to_input_directory}'@' ${output} 
#sed -i 's@FLE@'${input_file_name}'@' ${output} 
#sed -i 's@PRPO@'${path_github_repo}'@' ${output} 
#sed -i 's@x500@'${user_x500}'@' ${output} 
#sed -i 's@WK@'${desired_working_directory}'@' ${output} 
#sed -i 's@CRSMP@'${using_crossmap}'@' ${output} 
#sed -i 's@GNHRM@'${using_genome_harmonizer}'@' ${output}
#sed -i 's@KING@'${using_king}'@' ${output}
#sed -i 's@RPT@'${making_report}'@' ${output}
#sed -i 's@CSTQC@'${custom_qc}'@' ${output}
#sed -i 's@CSTANC@'${custom_ancestry}'@' ${output}
#sed -i 's@RFMX@'${using_rfmix}'@' ${output} 
#sed -i 's@LMDL@'${local_modules}'@' ${output}
