#!/bin/bash

work_dir=$(pwd)
echo ${work_dir}

cd ./src/QCReporter/

quarto render test_simpler.qmd
# quarto render simple.qmd #Currently not working
# quarto render QCReport.qmd #Currently not working

echo "Report has been successfully generated"
