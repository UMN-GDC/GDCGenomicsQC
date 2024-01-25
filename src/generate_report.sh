#!/bin/bash

cd ./src/QCReporter/

quarto render simple.qmd
# quarto render QCReport.qmd #Currently not working

echo "Report has been successfully generated"
