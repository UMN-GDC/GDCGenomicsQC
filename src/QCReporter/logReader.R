library(tidyverse)

#' @title Read the initial data from the log file'
#' @description Read the initial data from the log file'
#' @param filepath The path to the log file'
#' @return A list of named vectors characterizing the input and output data'
#' @export
#' @examples
#' extractLog("../QCtutorial/logs/wgas2.log")
extractLog<- function(filepath, plinkoption){
  file <- file(paste0(wd,filepath), "r")
  log <- readLines(file)
  close(file)


  test <- log[grep("people .* from .fam", log)]
  test1 <- str_extract_all(test, "[0-9]+")
  nSubjects <- as.numeric(test1[[1]][1])
  nMale <- as.numeric(test1[[1]][2])
  nFemale <- as.numeric(test1[[1]][3])

  # SNPsNSubjects()
  test <- log[grep("loaded from .bim", log)]
  nSNPs <- as.numeric(gsub("[^0-9]", "", test))

  initData <-c(nSubjects, nMale, nFemale, nSNPs)
  names(initData) <- c("InSubjects", "InMale", "InFemale", "InSNPs")


  # Passing QC
  # SNPs and subejcts
  test <- str_extract_all(log[grep("pass filters and QC", log)], "[0-9]+")
  nSNPs <- as.numeric(test[[1]][1])
  nSubjects <- as.numeric(test[[1]][2])
  
  #### Section for implementing switch to look for different things ####
  # All if string_to_find_a/string_to_find_b untested... 
  if(length(plinkoption) == 0) {plinkoption=5} 
  if (plinkoption == 1) { #--mind logs
    string_to_find_a="[0-9]+ people removed due to missing genotype data"
  }
  if (plinkoption ==2){ #--geno logs
    string_to_find_a="[0-9]+ variants removed due to missing genotype data"
  }

  if (plinkoption ==3) { #--check-sex logs
    string_to_find_a="[0-9]+ Xchr and [0-9]+ Ychr variant(s) scanned, [0-9]+ problems detected." 
  }
  
  if (plinkoption ==4) { #--maf logs
    string_to_find_a="[0-9]+ variants removed due to minor allele threshold(s)"
  } 
  if (plinkoption==5) {
    error_message= "Invalid plink option selected"
    stop(print(error_message))
  }
  
  string_to_find_b="[0-9]+ variants and [0-9]+ people pass filters and QC." #All logs have this in common!
  
  # case controls
  test2<- str_extract_all(log[grep('.* cases and .* are controls', log)], "[0-9]+")
  if(length(test2) == 0) {test2=NULL} #Extra patch... 
  if(is.null(test2) == FALSE) {
    nCases <- as.numeric(test2[[1]][1])
    nControls <- as.numeric(test2[[1]][2])
  } else {
    nCases = 000
    nControls = 000
  }
  outputData <- c(nSubjects, nCases, nControls, nSNPs)
  names(outputData) <- c("OutSubjects", "OutCases", "OutControls", "OutSNPs")
  # return nSubjects, nMale, nFemale, nSNPs
  return(c(initData, outputData))
}

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("A log to extract information from needs to be provided.", call.=FALSE)
} else if (length(args)>=1) {
  args[1] -> filename
  # default output file
  args[2] -> plink_option #mind, geno
  args[3] -> output_name
}
wd=getwd()
print(wd)
print(filename)
print(plink_option)
print(output_name)
final_name=paste0(wd, "/", output_name)
print(final_name)
plink_selected <- switch(plink_option,
       "mind"=1,       #=print("You have chosen to look for people missing genotype data logs."),
       "geno"=2,       #=print("You have chosen to look for variants missing genotype data logs."),
       "check-sex"=3,       #=print("You have chosen to look for comparisons between sex assingments logs."),
       "maf"=4)       #=print("You have chosen to look for SNPS with minor allele frequency logs."))

print(plink_selected)
table1=extractLog(filename, plink_selected)
table1
# write.table(table1, file = paste0(final_name))
#After extractLog need to save it as a table to be called within quarto


