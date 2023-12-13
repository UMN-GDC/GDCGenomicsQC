#!/bin/bash

Rscript logReader.R /sampleLogs/first_pass.log mind test1.txt #Works
Rscript logReader.R /sampleLogs/step2_temp.log geno test2.txt #Works
#Rscript logReader.R /sampleLogs/first_pass.log check-sex test3.txt
Rscript logReader.R /sampleLogs/SMILES_done_MAF.log maf test4.txt #Doesn't like the second part of it... 
#Rscript logReader.R /sampleLogs/gender_check.log mind test5.txt #This error message makes sense
#Rscript logReader.R /sampleLogs/bad_log.log mind test6.txt #This error message makes sense