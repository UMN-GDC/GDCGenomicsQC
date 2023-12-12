#!/bin/bash

Rscript logReader.R /sampleLogs/first_pass.log mind test1.txt
Rscript logReader.R /sampleLogs/first_pass.log geno test2.txt
Rscript logReader.R /sampleLogs/first_pass.log check-sex test3.txt
Rscript logReader.R /sampleLogs/first_pass.log maf test4.txt
Rscript logReader.R /sampleLogs/first_pass.log fakestuff test5.txt
