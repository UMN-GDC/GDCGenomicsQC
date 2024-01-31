#!/bin/bash

# Adding in the option for the user to choose what the file is called so PCA can be run independently
# and not always after QC steps are done

show_help() {
  echo "IMPORTANT!"
  echo "If you don't enter a filename using the --FILE flag this module assumes you ran main.sh previously"
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  --FILE <file_name>  Specify the file to perform PCA on. The file must be a .bed file"
  echo "  --help              Display this help message."
}

# Check for command line arguments
# if [ $# -eq 0 ]; then
#  echo "No arguments provided. Use --help for usage information."
#  exit 1
# fi
# Defaulting SWITCH
SWITCH=0

# Loop through the command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    --FILE)
        FILEUSING="$2"
        ((SWITCH++))
        echo "File chosen is ${FILEUSING%*}"
        shift 2 # Consume both the flag and its value
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


# Just in case they included the file extension
FILEUSING=${FILEUSING%*}
# echo "${SWITCH}"

# If they didn't include a file name
if [ ${SWITCH} == 0 ]; then
  FILE_LOOKUP=$(ls *_10.bed)  
  FILE_USING=${FILE_LOOKUP%.*}
  echo "${FILE_USING}"
fi

echo "(Step 10 of QC) Principle Component Analysis"

#We use fraposa to perform pca on the current data and a reference data 1000G with known population
#We superimpose the pcs of the current data onto the reference data to compare and predict the current data's population
git clone https://github.com/daviddaiweizhang/fraposa.git
mv ./fraposa/*.* ./ && rm -R fraposa
./commvar.sh 1000G ${FILE_USING} 1000G_comm ${FILE_USING}_11
./fraposa_runner.py --stu_filepref ${FILE_USING}_11 1000G_comm
./predstupopu.py 1000G_comm ${FILE_USING}_11
./plotpcs.py 1000G_comm ${FILE_USING}_11

# Removing junk files
rm -R __pycache__/