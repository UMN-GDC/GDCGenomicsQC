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

work_dir=$(pwd)
echo ${work_dir}

cd ./src/QCReporter/

final_location=${path_to_store_outputs}/results/${FILE}.pdf
path_to_data=${path_to_store_outputs}data/

# Script that alters where the data gets read from to be what gets provided in .qmd file #Alternatively could move all the files here temporarily...
# See the script I made for replacing lines of code in my TMS project
# Will need to replace line 37 with "path_to_data="/home/gdc/shared/GDC_pipeline/data/""
#### Script for altering a line of code here
file1=test_simpler.qmd
str1='path_to_data='\""${path_to_data}"'"'
# str1='path_to_data='\"''${path_to_data}''\"''
../replace_line.sh ${file1} 37 "${str1}"

quarto render test_simpler.qmd  

mv test_simpler.pdf ${final_location}



echo "Report has been successfully generated"
