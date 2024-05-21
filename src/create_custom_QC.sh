#!/bin/bash

# mkdir out
work_location=$1 #/home/miran045/shared/projects/Efield_modeling/experiments/msc
ids=$2 #sub-MSC04
full_path_to_simnibs_cifti_tools=$3
path_to_save_outputs=/scratch.global/miran045/shared/projects/Efield_modeling/gy_outputs/ #MSC




output=${work_location}/${i}_wrapper.sh
cp -v ${full_path_to_simnibs_cifti_tools}/wrapper_wip.sh ${output}
sed -i 's@G1@'${i}'@' ${output} #Default is 0.1
sed -i 's@M1@'${path_to_save_outputs}'@' ${output} #Default is 0.1
sed -i 's@G2@'${work_location}'/'${i}'/ses-01/list_subjs@' ${output} #Default is 0.1
sed -i 's@M2@'${work_location}'/'${i}'/ses-01/dlabel_tables/'${i}'_LRdlabel_table_native_no_mdwall.csv@' ${output} #Default is 0.1
sed -i 's@MAF1@'${full_path_to_simnibs_cifti_tools}'@' ${output} #Default is 0.1
sed -i 's@HWE1@'${var_hwe1}'@' ${output} #Default is 0.1
sed -i 's@HWE2@'${var_hwe2}'@' ${output} #Default is 0.1
