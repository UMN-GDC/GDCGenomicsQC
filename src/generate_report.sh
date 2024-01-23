#!/bin/bash

cd ./src/QCReporter/

quarto render QCReport.qmd

echo "Report has been successfully generated"
