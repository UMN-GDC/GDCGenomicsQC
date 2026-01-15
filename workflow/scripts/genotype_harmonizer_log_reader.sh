#!/bin/bash
workdir=$1

mkdir temp

# 1. Find all _aligned.log files
# 	* copy them to a temp folder
cd ${workdir}
cp *aligned.log temp/
cd ./temp/


# 2. Trim all _aligned.log files to be the last 16 rows
# Find all _aligned.log files in the current directory
aligned_files=$(find . -type f -name '*aligned.log')

# Iterate over each _aligned.log file
for file in ${aligned_files}; do
    # Get the filename without extension
    filename=$(basename "$file" "aligned.log")

    # Use the 'tail' command to get the last 16 rows of the file
    tail -n 16 "$file" > "${filename}_trimmed.log"

    echo "Trimmed $file to ${filename}_trimmed.log"
done

# Remove old _aligned.log files
rm *aligned.log


# 3. For each _trimmed.log Read the starting number of SNPs and which chromosome it is
# 4. Find out how many swapped SNPs 
# 	* Will need to sum this for convenience in reporting
# 5. Find the number of variants excluded during alignment phase

# Find all _trimmed.log files in the current directory
trimmed_files=$(find . -type f -name '*chr[0-9]*_trimmed.log')

# Initializing summary .txt file
echo Filename Chr Start_snps Swapped_snps Excluded_snps End_snps >genome_harmonizer_full_log.txt

# Iterate over each _trimmed.log file
for trimmed_file in $trimmed_files; do
    # Get the filename without extension
    filename=$(basename "$trimmed_file" "_trimmed.log")

    # Extract information from the file
    start_snps=$(cat ${trimmed_file} |  grep -oP 'Read \K\d+')
    chromosome=$(cat ${trimmed_file} | grep -oP 'chr\d+') 
    chromosome_b=$(echo "$chromosome" | tr -d 'chr')
    swapped_snps_1=$(cat ${trimmed_file} | grep -oP '(\d+,\d+)(?=\s+swapped)')
    swapped_snps_1b=$(echo "$swapped_snps_1" | tr -d ',')
    swapped_snps_2=$(cat ${trimmed_file} | grep -oP 'Swapped \K\d+')
    #summing swapped_snps
    total_swapped_snps=$((swapped_snps_1b + swapped_snps_2))
    excluded_variants=$(cat ${trimmed_file} | grep -oP 'total \K\d+,\d+')
    excluded_variantsb=$(echo "${excluded_variants}" | tr -d ',')
    end_snps=$(cat ${trimmed_file} | grep -oP 'Number of SNPs: \K\d+')
    excluded_variantsb=$((start_snps-end_snps))
    
    # Print the extracted information
    echo "File: ${filename}"
    echo "Chromosome: ${chromosome}"
    echo "Starting number of SNPs: ${start_snps}"
    echo "Number of swapped SNPs: ${total_swapped_snps}"
    echo "Number of variants excluded during alignment: ${excluded_variantsb}"
    echo "End number of SNPs: ${end_snps}"
    echo "-------------------------"
    
    # Save the information in a .txt file
    echo ${filename} ${chromosome_b} ${start_snps} ${total_swapped_snps} ${excluded_variantsb} ${end_snps} >>genome_harmonizer_full_log.txt    
done

cat genome_harmonizer_full_log.txt
# Moving desired output into prior folder
mv genome_harmonizer_full_log.txt ..

# Removing the temporary files from this step
cd ..
rm -R temp/


