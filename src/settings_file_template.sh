#!/bin/bash

/home/gdc/shared/GDC_pipeline/GDCGenomicsQC/src/settings_file_reader.sh \
--path_to_input_directory /path/to/study_here \
--input_file_name sample_study \
--path_github_repo /home/gdc/shared/GDC_pipeline/GDCGenomicsQC \
--user_x500 samp300 \
--desired_working_directory /path_to/working_area/stores_outputs/here \
--using_crossmap 1 \
--using_genome_harmonizer 1 \
--making_report 1 \
--custom_qc 0

