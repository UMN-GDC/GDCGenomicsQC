#!/bin/bash

Rscript logReader.R /sampleLogs/first_pass.log #Updated error detecting works
Rscript logReader.R /sampleLogs/first_pass.log mind #Default output file name works
Rscript logReader.R /sampleLogs/first_pass.log mind test1.txt #Works
Rscript logReader.R /sampleLogs/step2_temp.log geno test2.txt #Works
Rscript logReader.R /sampleLogs/gender_check.log check-sex test3.txt #Works
Rscript logReader.R /sampleLogs/SMILES_done_MAF.log maf test4.txt #Works
Rscript logReader.R /sampleLogs/Step7_temp1.log filter-founders test5.txt #Works 
Rscript logReader.R /sampleLogs/Step5_temp.log hwe test6.txt #Works
Rscript logReader.R /sampleLogs/gender_check.log mind test5.txt #For testing error messages.. This error message makes sense 
Rscript logReader.R /sampleLogs/bad_log.log mind test6.txt #For testing error messages.. This error message makes sense

#Plink options I am unsure of... Probably not used in the pipeline...
# --mpheno 4 [SMILES_pheno_updated.log
# same with [SMILES_matched_to_b38_forward.log]
# [SMILES_comm.log] <-- FRAPOSA
# --het [pruned_data.log] <-- doesn't really give new information
# --extract [pihat_min0.2_in_founders.log] <-- unsure what it adds
# indepSNP.log <-- This one looks useful just difficult to incorporate. Do we want each chromosome information?
# Unsure about [hgdp_v1.log]
# 