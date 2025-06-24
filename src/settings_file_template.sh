#!/bin/bash
path_github_repo=Put/path/GDCGenomicsQC/github/repo

${path_github_repo}/src/settings_file_reader.sh \
--path_to_input_directory /path/to/study_here \
--input_file_name sample_study \
--path_github_repo ${path_github_repo} \
--user_x500 samp300 \
--desired_working_directory /path_to/working_area/stores_outputs/here \
--using_crossmap 1 \
--using_genome_harmonizer 1 \
--making_report 1 \
--custom_qc 0 \
--using_rfmix 1 \
--use_primus 0

