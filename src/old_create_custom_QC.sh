#!/bin/bash


initial_marker_filtering=$1 
initial_sample_filtering=$2
ultimate_marker_filtering=$3
ultimate_sample_filtering=$4
maf_filtering=$5
hwe_controls=$6
hwe_cases=$7
work_location=$8
full_path_to_repo=$9


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
