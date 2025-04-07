#!/bin/bash

show_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  --FILE <file_name>  Specify the file to gather QC information after running main.sh. Must be a .bed file"
  echo "  --PATHTOSTOREOUTPUTS  Specify the full path to where you would like the outputs of this pipeline to go"
  echo "  --help              Display this help message."
}

# Check for command line arguments
if [ $# -eq 0 ]; then
  echo "No arguments provided. Use --help for usage information."
  exit 1
fi

# Loop through the command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    --FILE)
        FILE="$2"
        echo "File chosen is $FILE"
        shift 2 # Consume both the flag and its value
      ;;
    --PATHTOSTOREOUTPUTS)
        path_to_store_outputs="$2"
        echo "Path provided is $path_to_store_outputs"
        shift 2
      ;;   
    --help)
        show_help
        exit 0
      ;;
    *)
      echo "Unrecognized option: $key"
      show_help
      exit 1
      ;;
  esac
done


# Just in case the file extension is included
FILE=${FILE%*}

current_dir=$(pwd)
echo "You are currently in: $current_dir"

place_to_store_data=${path_to_store_outputs}
mkdir -p ${place_to_store_data}

path_to_logReader=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/QCReporter/logReader.R
path_to_logReader_extended=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/QCReporter/logReader_extended.R

#Gathers all the information from logs and puts them into tables for later use
Rscript ${path_to_logReader} ${path_to_store_outputs}/${FILE}1.log geno QC2_geno.txt ${place_to_store_data}
Rscript ${path_to_logReader} ${path_to_store_outputs}/${FILE}2.log mind QC3_mind.txt ${place_to_store_data}
Rscript ${path_to_logReader} ${path_to_store_outputs}/${FILE}3.log geno QC4_geno.txt ${place_to_store_data}
Rscript ${path_to_logReader} ${path_to_store_outputs}/${FILE}4.log mind QC5_mind.txt ${place_to_store_data}
#Rscript ${path_to_logReader} ${path_to_store_outputs}/gender_check.log check-sex QC6_sex_check.txt ${place_to_store_data}
Rscript ${path_to_logReader} ${path_to_store_outputs}/${FILE}6.log maf QC_7_maf.txt ${place_to_store_data} 
Rscript ${path_to_logReader} ${path_to_store_outputs}/${FILE}7a.log hwe QC_8_hwe.txt ${place_to_store_data}
Rscript ${path_to_logReader} ${path_to_store_outputs}/${FILE}7.log hwe QC_8b_hwe.txt ${place_to_store_data}
#Rscript ${path_to_logReader} ${path_to_store_outputs}/${FILE}9a.log filter-founders QC_9_filter-founders.txt ${place_to_store_data}

Rscript ${path_to_logReader_extended} ${path_to_store_outputs}/indepSNP.log indep-pairwise QC_indep_pairwise.txt ${place_to_store_data} 

