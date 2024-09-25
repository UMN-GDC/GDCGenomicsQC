#!/bin/bash


path_to_input_directory=$1 
input_file_name=$2
path_github_repo=$3
user_x500=$4
desired_working_directory=$5
using_crossmap=$6
using_genome_harmonizer=$7
using_rfmix=$8
making_report=$9
custom_qc=${10}

output=${desired_working_directory}/${input_file_name}_wrapper.sh
mkdir -p ${desired_working_directory}
cp -v ${path_github_repo}/src/main_template.sh ${output}
sed -i 's@PND@'${path_to_input_directory}'@' ${output} 
sed -i 's@FLE@'${input_file_name}'@' ${output} 
sed -i 's@PRPO@'${path_github_repo}'@' ${output} 
sed -i 's@x500@'${user_x500}'@' ${output} 
sed -i 's@WK@'${desired_working_directory}'@' ${output} 
sed -i 's@CRSMP@'${using_crossmap}'@' ${output} 
sed -i 's@GNHRM@'${using_genome_harmonizer}'@' ${output}
sed -i 's@RFMX@'${using_rfmix}'@' ${output} 
sed -i 's@RPT@'${making_report}'@' ${output}
sed -i 's@CSTQC@'${custom_qc}'@' ${output}
