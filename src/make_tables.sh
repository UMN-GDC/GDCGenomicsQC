#!/bin/bash

show_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  --FILE <file_name>  Specify the file to process. Must be a .bed file"
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
    --help)
        show_help
        exit 0
      ;;
          *) 
            echo "Please enter a valid option (Yes/No) when choosing to include the UPDTSO flag."
            exit 0
          ;;
        esac
        shift 2
      ;;
    *)
      echo "Unrecognized option: $key"
      show_help
      exit 1
      ;;
  esac
done


# Just in case they included the file extension
FILE=${FILE%*}


#Gathers all the information from logs and puts them into tables for later use
Rscript logReader.R ./${FILE}_1.log geno QC2_geno.txt
Rscript logReader.R ./${FILE}_2.log mind QC3_mind.txt
Rscript logReader.R ./${FILE}_3.log geno QC4_geno.txt
Rscript logReader.R ./${FILE}_4.log mind QC5_mind.txt
Rscript logReader.R ./gender_check.log check-sex QC6_sex_check.txt
Rscript logReader.R ./${FILE}_6_MAF.log maf QC_7_maf.txt
## Made it to QC step 8 ##


#for troubleshooting and reference
Rscript logReader.R /sampleLogs/first_pass.log mind #Default output file name works
Rscript logReader.R /sampleLogs/first_pass.log mind test1.txt #Works
Rscript logReader.R /sampleLogs/step2_temp.log geno test2.txt #Works
Rscript logReader.R /sampleLogs/gender_check.log check-sex test3.txt #Works
Rscript logReader.R /sampleLogs/SMILES_done_MAF.log maf test4.txt #Works
Rscript logReader.R /sampleLogs/Step7_temp1.log filter-founders test5.txt #Works 
Rscript logReader.R /sampleLogs/Step5_temp.log hwe test6.txt #Works
