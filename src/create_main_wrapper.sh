#!/bin/bash


path_to_input_directory=$1 
input_file_name=$2
path_github_repo=$3
user_x500=$4
desired_working_directory=$5
using_crossmap=$6
using_genome_harmonizer=$7
using_king=$8
using_rfmix=$9
making_report=${10}
custom_qc=${11}
combine_related=${12}
custom_ancestry=${13}
CHECK_SEX=${14}

# Check environment if we are currently running inside Apptainer
if [[ -n "$APPTAINER_NAME" || -n "$SINGULARITY_NAME" ]]; then
    echo "Environment: Inside Container ($APPTAINER_NAME)"
    path_github_repo=/app/GDCGenomicsQC
    template=${path_github_repo}/src/main_template_image.sh
else
    echo "Environment: Native Host"
    template=${path_github_repo}/src/main_template.sh

fi

mkdir -p ${desired_working_directory}
cp -v ${template} ${output}
sed -i 's@PND@'${path_to_input_directory}'@' ${output} 
sed -i 's@FLE@'${input_file_name}'@' ${output} 
sed -i 's@PRPO@'${path_github_repo}'@' ${output} 
sed -i 's@x500@'${user_x500}'@' ${output} 
sed -i 's@WK@'${desired_working_directory}'@' ${output} 
sed -i 's@CRSMP@'${using_crossmap}'@' ${output} 
sed -i 's@GNHRM@'${using_genome_harmonizer}'@' ${output}
sed -i 's@KING@'${using_king}'@' ${output}
sed -i 's@RPT@'${making_report}'@' ${output}
sed -i 's@CSTQC@'${custom_qc}'@' ${output}
sed -i 's@COMB@'${combine_related}'@' ${output}
sed -i 's@CSTANC@'${custom_ancestry}'@' ${output}
sed -i 's@RFMX@'${using_rfmix}'@' ${output} 
sed -i 's@CHECK_SEXX@'${CHECK_SEX}'@' ${output} 
