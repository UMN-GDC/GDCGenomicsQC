#!/bin/bash

show_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "--path_to_input_directory		Provide the path to where the bim/bed/fam data is stored"
  echo "--input_file_name 				Provide the common name that ties the bim/bed/fam files together"
  echo "--path_github_repo 				Provide the path to the GDCGenomicsQC pipeline repository"
  echo "--user_x500						Provide your x500 samp213@umn.edu so you may receive email updates regarding sbatch submissions"
  echo "--desired_working_directory		Provide a path to where you'd like the outputs to be stored"
  echo "--using_crossmap				Enter '1' for if you would like to update your reference genome build from GRCh37 to GRCh38"
  echo "--using_genome_harmonizer 		Enter '1' if you would like to update strand allignment by using genome harmonizer"
  echo "--making_report					Enter '1' if you would like an automated report to be generated of the qc steps and what was changed"
  echo "--custom_qc						Enter '1' if you would like to use your own settings for the qc steps such as marker and sample filtering"
  echo "								When providing this flag you will need to answer all of the questions prompted by the terminal"
  echo "Default settings: The pipeline by default if flags are not provided will use crossmap, genome harmonizer and will generate the automated reports"
  echo "  --help              Display this help message."
}


# Check for command line arguments
if [ $# -eq 0 ]; then
  echo "No arguments provided. Use --help for usage information."
  exit 1
fi

#Default settings
NYS=1
output_folder=$(pwd)
path_github_repo=/home/gdc/shared/GDC_pipeline/GDCGenomicsQC
using_crossmap=1
using_genome_harmonizer=1
making_report=1
custom_qc=0
desired_working_directory=/scratch.global/gdc

# Loop through the command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --path_to_input_directory)
        path_to_input_directory="$2"
        echo "Provided path to the input directory is ${path_to_input_directory}"
        shift 2 # Consume both the flag and its value
      ;;
    --input_file_name)
        input_file_name="$2"
        echo "Plink bed formatted file name is ${input_file_name}"
        shift 2 # Consume both the flag and its value
      ;;
    --path_github_repo)
        path_github_repo="$2"
        # Check if the user provided a valid path to the repo
        if [ ! -f "${path_github_repo}/README.md" ]; then
            echo "Incorrect directory provided: ${path_github_repo}"
            exit 1
        fi
        shift 2
      ;;
    --user_x500)
        user_x500="$2"
        echo "Your umn id is ${user_x500}"
        shift 2 # Consume both the flag and its value
      ;;
    --desired_working_directory)
        desired_working_directory="$2"
        echo "Your desired working directory for scripts and outputs is ${desired_working_directory}"
        shift 2 
      ;;
    --using_crossmap)
        using_crossmap="$2"
        if [ ${using_crossmap} -eq 1 ]; then
          echo "You have chosen to update your reference genome build from GRCh37 to GRCh38"
        fi
        shift 2
      ;;
    --using_genome_harmonizer)
        using_genome_harmonizer="$2"
        if [ ${using_genome_harmonizer} -eq 1 ]; then
          echo "You have chosen to update your strand allignment by using genome harmonizer"
        fi
        shift 2 
      ;;
    --making_report)
        making_report="$2"
        if [ ${making_report} -eq 1 ]; then
          echo "You have chosen to utilize the automated report feature"
        fi
        shift 2
      ;;
    --custom_qc)
        custom_qc="$2"
        if [ ${custom_qc} -eq 1 ]; then
          echo "You have chosen to customize the standard qc steps, please answer all of the following questions"
          source ${path_github_repo}/src/gather_custom_qc_inputs.sh ${path_github_repo} ${desired_working_directory}
        fi
        shift 2
      ;;
    *)
      echo "Unrecognized option: $key"
      show_help
      exit 1
      ;;
  esac
done


## Making personalized wrapper for submitting GenomicsQCPipeline in users specified working_directory
source ${path_github_repo}/src/create_main_wrapper.sh \
${path_to_input_directory} \
${input_file_name} \
${path_github_repo} \
${user_x500} \
${desired_working_directory} \
${using_crossmap} \
${using_genome_harmonizer} \
${making_report} \
${custom_qc}


