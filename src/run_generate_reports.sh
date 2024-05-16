#!/bin/bash

module load plink
module load python3/3.9.3_anaconda2021.11_mamba
module load R

# Could add in a path_to_repo argument so that the paths to functions won't be hardcoded anymore
working_directory=$1 #Only argument needed for this process to run
# working_directory=/home/gdc/shared/GDC_pipeline/results/Needed_files_for_report/SMILES_GDA


array_location_base=($(ls -d ${working_directory}/*/))
num_elements=${#array_location_base[@]}

# array_location=(Full EUR AMR AFR SAS)
# filepreffix=(mixed.ethnic.QC EUR.QC AMR.QC AFR.QC SAS.QC)
path_to_replace_line_function=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/replace_line.sh

path_to_qmd=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/QCReporter
path_to_gen_all_reports=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/QCReporter
file1=${path_to_qmd}/updated_report_template.qmd #Full path to report
file2=${path_to_qmd}/ancestry_report_template.qmd

for ((i=0; i<${num_elements}; i++)); do
    array_location[i]=${array_location_base[$i]%/} #This gets the full paths to each directory
    filepreffix_test[i]=${array_location[$i]##*/}
    filepreffix[i]=${filepreffix_test[i]}.QC
    echo ${filepreffix[i]}
    echo ${filepreffix_test[i]}
    
    path_to_store_outputs=${array_location[i]}
    ${path_to_gen_all_reports}/generate_all_reports.sh --FILE ${filepreffix[i]} --PATHTOSTOREOUTPUTS ${path_to_store_outputs} 

    final_location=${path_to_store_outputs}/results/${filepreffix[i]}.pdf
    mkdir -p ${path_to_store_outputs}/results

## Changes where the qmd looks for the data... 
    gender_file_name=$(ls ${path_to_store_outputs}/*.sexcheck) #Returns the full path
    pushd ${path_to_store_outputs}
    gender_file_name=$(ls *.sexcheck)
    popd
    path_read_files=${path_to_store_outputs}/
    
    ## Making QC report based on saved template
    output_qmd_1=${path_to_store_outputs}/results/${filepreffix[i]}.qmd
    cp -v ${file1} ${output_qmd_1}
    sed -i 's@PATH@'${path_read_files}'@' ${output_qmd_1}
    sed -i 's@NAME@'${gender_file_name}'@' ${output_qmd_1}

    ## Generating report
    quarto render ${output_qmd_1} 

    
    if [ -f "${final_location}" ]; then
      rm ${output_qmd_1}
      echo "Report has been successfully generated for ${array_location[i]}"
      echo "Removing temp quarto document"
    else
      echo "Skipping removal of temp file since an error occured while generating report"
    fi
    
    
    if [ "${filepreffix_test[i]}" == "full" ]; then
      ## For full ancestry directory
      pushd ${path_to_store_outputs}
      fraposa_log_file=$(ls *.unrelated.comm.popu)
      popd
      frapose_png=$(ls ${path_to_store_outputs}/*.unrelated.comm*.png)
      
      output_qmd_2=${path_to_store_outputs}/results/ancestry_report.qmd
      cp -v ${file2} ${output_qmd_2}
      sed -i 's@PATH@'${path_read_files}'@' ${output_qmd_2}
      sed -i 's@NAME@'${fraposa_log_file}'@' ${output_qmd_2}
      sed -i 's@SED@'${frapose_png}'@' ${output_qmd_2} 
    

      quarto render ${output_qmd_2} 
      final_location_2=${path_to_store_outputs}/results/ancestry_report.pdf

      if [ -f "${final_location_2}" ]; then
        rm ${output_qmd_2}
        echo "Ancestry report has been successfully generated for ${array_location[i]}"
        echo "Removing temp quarto document"
      else
        echo "Skipping removal of temp file since an error occured while generating report"
      fi
    fi

done
