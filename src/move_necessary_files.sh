#!/bin/bash


show_help() {
  echo "IMPORTANT!"
  echo "If you don't enter a filename using the --FILE flag this module assumes you ran main.sh previously"
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  --FILE <file_name>  Specify the text file which has a list of files to remove from the repo after each run."
  echo "  --help              Display this help message."
}

# Check for command line arguments
# if [ $# -eq 0 ]; then
#  echo "No arguments provided. Use --help for usage information."
#  exit 1
# fi
# Defaulting SWITCH
SWITCH=0
switch=0

# Loop through the command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    --FILE)
        text_file="$2"
        ((SWITCH++))
        echo "File chosen is ${text_file}"
        shift 2 # Consume both the flag and its value
      ;;
    --END_LOCATION)
        final_location="$2"
        ((switch++))
        echo "Spot to move files to is ${final_location}"
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


# If they didn't include a file name
## Specify the text file containing the list of files to be removed

if [ ${SWITCH} == 0 ]; then
  text_file=./data/default_list_of_files_to_move.txt
fi

if [ ${switch} == 0 ]; then
  final_location=/home/gdc/shared/GDC_pipeline/logs/
fi

# Check if the text file exists
if [ -f "${text_file}" ]; then
    # Loop through each line in the text file and remove the corresponding file
    while IFS= read -r file_to_move; do
        if [ -e "${file_to_move}" ]; then
            mv "${file_to_move}" ${final_location}/
            echo "Moved: ${file_to_move}"
        else
            echo "File not found: ${file_to_move}"
        fi
    done < "${text_file}"
else
    echo "Error: Text file '${text_file}' not found."
fi

echo "All done cleaning up the repository!"
