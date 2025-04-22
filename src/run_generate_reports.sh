#!/bin/bash

# Usage: ./this_script.sh <working_directory> <path_to_repo>

set -uo pipefail

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <working_directory> <path_to_repo>"
    exit 1
fi

working_directory="$1"
path_to_repo="$2"

path_to_qmd="${path_to_repo}/src/QCReporter"
template_main="${path_to_qmd}/updated_report_template.qmd"
template_ancestry="${path_to_qmd}/ancestry_report_template.qmd"
generate_script="${path_to_repo}/src/QCReporter/generate_all_reports.sh"

generate_report() {
    local subject_dir="$1"
    local target_dir="$2"
    local path_to_repo="$3"
    local file_prefix="${target_dir}.QC"
    local output_dir="${subject_dir}/results"
    local output_qmd="${output_dir}/${file_prefix}.qmd"
    local final_pdf="${output_dir}/${file_prefix}.pdf"

    mkdir -p "$output_dir"

    # Run generate_all_reports
    "$generate_script" --FILE "$file_prefix" --PATHTOSTOREOUTPUTS "$subject_dir" --path_to_repo "${path_to_repo}"

    # Find gender file
    local gender_file_path
    gender_file_path=$(find "$subject_dir" -maxdepth 1 -name "*.sexcheck" | head -n 1)
    if [ -z "$gender_file_path" ]; then
        echo "No .sexcheck file found in $subject_dir. Skipping."
        return
    fi
    local gender_file_name
    gender_file_name=$(basename "$gender_file_path")

    # Prepare the report file
    cp -v "$template_main" "$output_qmd"
    sed -i -e "s@PATH@${subject_dir}/@" \
           -e "s@NAME@${gender_file_name}@" "$output_qmd"

    # Render the report
    if quarto render "$output_qmd"; then
        rm -v "$output_qmd"
        echo "Report successfully generated for $target_dir"
    else
        echo "Error generating report for $target_dir"
    fi
}

generate_ancestry_report() {
    local subject_dir="$1"
    local output_dir="${subject_dir}/results"
    local output_qmd="${output_dir}/ancestry_report.qmd"
    local final_pdf="${output_dir}/ancestry_report.pdf"

    mkdir -p "$output_dir"

    local fraposa_log
    fraposa_log=$(find "$subject_dir" -maxdepth 1 -name "*.unrelated.comm.popu" | head -n 1)
    local pc1_pc2
    pc1_pc2=$(find "$subject_dir" -name "PC1*PC2.png" | head -n 1)
    local pc1_pc3
    pc1_pc3=$(find "$subject_dir" -name "PC1*PC3.png" | head -n 1)
    local pc2_pc3
    pc2_pc3=$(find "$subject_dir" -name "PC2*PC3.png" | head -n 1)

    if [[ -z "$fraposa_log" || -z "$pc1_pc2" || -z "$pc1_pc3" || -z "$pc2_pc3" ]]; then
        echo "L Missing ancestry components in $subject_dir. Skipping ancestry report."
        return
    fi

    cp -v "$template_ancestry" "$output_qmd"
    sed -i -e "s@PATH@${subject_dir}/@" \
           -e "s@NAME@$(basename "$fraposa_log")@" \
           -e "s@PL1P@${pc1_pc2}@" \
           -e "s@PL2P@${pc1_pc3}@" \
           -e "s@PL3P@${pc2_pc3}@" "$output_qmd"

    if quarto render "$output_qmd"; then
        rm -v "$output_qmd"
        echo "Ancestry report successfully generated for $(basename "$subject_dir")"
    else
        echo "Error generating ancestry report for $(basename "$subject_dir")"
    fi
}

# Main loop
for dir in "$working_directory"/*/; do
    dir="${dir%/}"
    target_dir="${dir##*/}"

    [[ "$target_dir" == "temp" ]] && continue

    echo "Processing subject: $target_dir"

    generate_report "$dir" "$target_dir" "$path_to_repo"

    if [[ "$target_dir" == "full" ]]; then
        generate_ancestry_report "$dir"
    fi
done