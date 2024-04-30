#!/bin/bash

module load plink
module load python3/3.9.3_anaconda2021.11_mamba
module load R

working_directory=/home/gdc/shared/GDC_pipeline/results/Needed_files_for_report/SMILES_GSA

array_location_base=($(ls -d ${working_directory}/*/))
num_elements=${#array_location_base[@]}



for ((i=0;i<${num_elements}; i++)); do
    (array_location[i]=${array_location_base[$i]%/} #This gets the full paths to each directory
    filepreffix_test[i]=${array_location[$i]##*/}
    filepreffix[i]=${filepreffix_test[i]}.QC
    echo ${filepreffix[i]}) &
done
wait

# array_location=(Full EUR AMR AFR SAS)
# filepreffix=(mixed.ethnic.QC EUR.QC AMR.QC AFR.QC SAS.QC)
path_to_replace_line_function=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/replace_line.sh

path_to_qmd=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/QCReporter
path_to_gen_all_reports=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/QCReporter
file1=${path_to_qmd}/updated_report.qmd #Full path to report

for ((i=0; i<${num_elements}; i++)); do
    path_to_store_outputs=${array_location[i]}
    ${path_to_gen_all_reports}/generate_all_reports.sh --FILE ${filepreffix[i]} --PATHTOSTOREOUTPUTS ${path_to_store_outputs} 

    final_location=${path_to_store_outputs}/results/${filepreffix[i]}.pdf
    mkdir -p ${path_to_store_outputs}/results
    # path_to_data=${path_to_store_outputs}

## Changes where the qmd looks for the data... 
  # Need to make this updated_report.qmd file from the current working test_simpler.qmd document##
    gender_file_name=$(ls ${path_to_store_outputs}/*.sexcheck) #Still returns the full path
    pushd ${path_to_store_outputs}
    gender_file_name=$(ls *.sexcheck)
    popd
    path_read_files=${path_to_store_outputs}/
    str1='path_to_data='\""${path_read_files}"'"'
    str2='gender_file_name='\""${gender_file_name}"'"'
    # str1='path_to_data='\"''${path_to_data}''\"''
    ${path_to_replace_line_function} ${file1} 39 "${str1}"
    ${path_to_replace_line_function} ${file1} 40 "${str2}"

    quarto render ${file1} -P path_to_data:${path_to_data} -P gender_file_name:'plink.sexcheck' #The parameters here don't really work yet

    mv -v ${path_to_qmd}/updated_report.pdf ${final_location}

    echo "Report has been successfully generated"

done
