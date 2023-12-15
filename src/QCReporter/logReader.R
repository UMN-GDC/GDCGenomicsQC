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
  
  num_args=1
  #### Section for implementing switch to look for different things ####
  # Testing if this works for plinkoption ==1 
  if(length(plinkoption) == 0) {plinkoption=7} 
  if (plinkoption == 1) { #--mind logs
    string_to_find_a="[0-9]+ people removed due to missing genotype data"
  }
  if (plinkoption ==2){ #--geno logs
    string_to_find_a="[0-9]+ variants removed due to missing genotype data"
  }

  if (plinkoption ==3) { #--check-sex logs
    string_to_find_a=".* Xchr and .* Ychr variant.* scanned, .* problems detected." 
    num_args=3
  }
  
  if (plinkoption ==4) { #--maf logs 
    string_to_find_a=".*variants removed due to minor allele threshold.*"
  } 
  if (plinkoption ==5) { #--filter-founders logs 
    string_to_find_a=".* people removed due to founder status .*"
    num_args=2
  }
  if (plinkoption ==6) { #--hwe logs 
    string_to_find_a=".*variants removed due to Hardy-Weinberg exact test."
  }
  if (plinkoption==7) { #Default
    error_message= "Invalid plink option selected"
    stop(print(error_message))
  }
  
  string_to_find_b="[0-9]+ variants and [0-9]+ people pass filters and QC." #All logs have this in common!
  string_to_find_c="Before main variant filters, .* founders and .* nonfounders present."
  
  # Extracting what was changed
  ## Works for option ==1, 2, 4
  ## Now for option == 3

  test2<- str_extract_all(log[grep(string_to_find_a, log)], "[0-9]+")
  if(length(test2) == 0) {
    error_message= "Invalid log file provided for this plink option selected"
    stop(print(error_message))
  } else{
    nRemoved <- as.numeric(test2[[1]][1])
  }
  if(num_args ==3) {
    nY <- as.numeric(test2[[1]][2])
    nProblems <- as.numeric(test2[[1]][3])
  }
  if(num_args ==2) {
    test4 <- str_extract_all(log[grep(string_to_find_c, log)], "[0-9]+")
    if(length(test4) == 0) {
      error_message= "Invalid log file provided for this plink option selected"
      stop(print(error_message))
    }
    nFounder <- as.numeric(test4[[1]][1])
    nNoFounder <- as.numeric(test4[[1]][2])
  }
        

 
  test3<- str_extract_all(log[grep(string_to_find_b, log)], "[0-9]+")
  if(length(test3)==0){
    error_message_2="Bad string_to_find_b"
    stop(print(error_message_2))
  }
  nSNPS2 <- as.numeric(test3[[1]][1])
  nSubjs2 <- as.numeric(test3[[1]][2])
  outputData <- c(nSubjs2, nRemoved, nSNPS2)
  ## Would be a good idea to have the column names for outputData be different for each plink option (1,2,4)
  names(outputData) <- c("OutSubjects", "NumRemoved", "OutSNPs")
  
  if(num_args ==3) { #Output for plink option 3
  outputData <- c(nSubjs2, nRemoved, nY, nProblems, nSNPS2)
  names(outputData) <- c("OutSubjects", "NumX", "NumY", "NumProblems", "OutSNPs")
  }
  
  if(num_args ==2) { #Output for plink option 5
  outputData <- c(nSubjs2, nRemoved, nFounder, nNoFounder, nSNPS2)
  names(outputData) <- c("OutSubjects", "NumRemoved", "NumFounder", "NumNonFounder", "OutSNPs")
  }
  
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
       "maf"=4,       #=print("You have chosen to look for SNPS with minor allele frequency logs."))
       "filter-founders" =5,
       "hwe" = 6)

print(plink_selected)
table1=extractLog(filename, plink_selected)
table1
# write.table(table1, file = paste0(final_name))
#After extractLog need to save it as a table to be called within quarto


