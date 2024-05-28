#!/bin/bash
read -p "Please provide the path to the GDCGenomicsQC repository: " full_path_to_repo

echo "The default value for initial marker filtering is 0.1"
read -p "Enter your desired initial marker filtering threshold: " initial_marker_filtering

echo "The default value for initial sample filtering is 0.1" 
read -p "Enter your desired initial sample filtering threshold: " initial_sample_filtering

echo "The default value for ultimate marker filtering is 0.02"
read -p "Enter your desired ultimate marker filtering threshold: " ultimate_marker_filtering

echo "The default value for ultimate sample filtering is 0.02"
read -p "Enter your desired ultimate sample filtering threshold: " ultimate_sample_filtering

echo "The default value for minor allele frequency filtering is 0.01"
read -p "Enter your desired minor allele frequency filtering threshold: " maf_filtering

echo "By default markers where the hwe p-values are less than 1e-10 are removed for cases"
read -p "Enter your desired lower limit for hwe p-values for cases: " hwe_controls

echo "For controls markers are removed with hwe p-values with less than 1e-6"
read -p "Enter your desired lower limit for hwe p-values for controls " hwe_cases

read -p "Please provide where you would like your custom_qc.SLURM script to be stored: " work_location


output=${work_location}/custom_qc.SLURM
mkdir -p ${work_location}
cp -v ${full_path_to_repo}/src/QC_template.SLURM ${output}
sed -i 's@G1@'${initial_marker_filtering}'@' ${output} #Default is 0.1
sed -i 's@M1@'${initial_sample_filtering}'@' ${output} #Default is 0.1
sed -i 's@G2@'${ultimate_marker_filtering}'@' ${output} #Default is 0.02
sed -i 's@M2@'${ultimate_sample_filtering}'@' ${output} #Default is 0.02
sed -i 's@MAF1@'${maf_filtering}'@' ${output} #Default is 0.01
sed -i 's@HWE1@'${hwe_controls}'@' ${output} #Default is 1e-6
sed -i 's@HWE2@'${hwe_cases}'@' ${output} #Default is 1e-10
