#' @title Read the initial data from the log file'
#' @description Read the initial data from the log file'
#' @param filepath The path to the log file'
#' @return A list of named vectors characterizing the input and output data'
#' @export
#' @examples
#' extractLog("../QCtutorial/logs/wgas2.log")

suppressMessages(library(tidyverse)) 

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
  
  #Initializing indicators
  out_label="NotUsed"
  num_args=1
  
  # Section for implementing switch to look for different things  
  if(length(plinkoption) == 0) {plinkoption=7} 
  if (plinkoption == 1) { #--mind logs
    string_to_find_a="[0-9]+ people removed due to missing genotype data"
    out_label="NumPeopleRemoved"
  }
  if (plinkoption ==2){ #--geno logs
    string_to_find_a="[0-9]+ variants removed due to missing genotype data"
    out_label="NumVariantsRemoved"
  }
  if (plinkoption ==3) { #--check-sex logs
    string_to_find_a=".* Xchr and .* Ychr variant.* scanned, .* problems detected." 
    num_args=3
  }
  if (plinkoption ==4) { #--maf logs 
    string_to_find_a=".*variants removed due to minor allele threshold.*"
    out_label="NumVariantsRemoved"
  } 
  if (plinkoption ==5) { #--filter-founders logs 
    string_to_find_a=".* people removed due to founder status .*"
    num_args=2
  }
  if (plinkoption ==6) { #--hwe logs 
    string_to_find_a=".*variants removed due to Hardy-Weinberg exact test."
    out_label="NumVariantsRemoved"
  }
  if (plinkoption==7) { #Default
    error_message= "Invalid plink option selected"
    stop(print(error_message))
  }
  
  string_to_find_b="[0-9]+ variants and [0-9]+ people pass filters and QC." #All logs have this in common!
  # for use in plinkoption 5
  string_to_find_c="Before main variant filters, .* founders and .* nonfounders present."
  
  # Extracting what was changed
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
        

  # Retreiving surviving #SNPs and participants
  test3<- str_extract_all(log[grep(string_to_find_b, log)], "[0-9]+")
  if(length(test3)==0){
    error_message_2="Bad string_to_find_b"
    stop(print(error_message_2))
  }
  nSNPS2 <- as.numeric(test3[[1]][1])
  nSubjs2 <- as.numeric(test3[[1]][2])
  outputData <- c(nSubjs2, nRemoved, nSNPS2)
  
  # Allows for different column names between plink options
  names(outputData) <- c("OutSubjects", out_label, "OutSNPs")
  
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

# test if there is at least two arguments: if not, return an error
if (length(args)<=1) {
  stop("A log to extract information from and which plink option need to be provided.", call.=FALSE)
} else if (length(args)>=2) {
  args[1] -> filename
  # default output file
  args[2] -> plink_option #mind, geno
  if(length(args)==2){
    output_name <- "extractLog_output.txt"
  } else {
  args[3] -> output_name
  }
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
print(" ")
write.table(table1, file = paste0(final_name), quote = F, row.names = T, col.names = F)
#After extractLog need to save it as a table to be called within quarto


=======
library(tidyverse)

#' @title Read the initial data from the log file'
#' @description Read the initial data from the log file'
#' @param filepath The path to the log file'
#' @return A list of named vectors characterizing the input and output data'
#' @export
#' @examples
#' extractLog("../QCtutorial/logs/wgas2.log")
extractLog<- function(filepath){
  file <- file(filepath, "r")
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
  
  
  # case control
  test2<- str_extract_all(log[grep(".* cases and .* are controls", log)], "[0-9]+")
  nCases <- as.numeric(test2[[1]][1])
  nControls <- as.numeric(test2[[1]][2])

  outputData <- c(nSubjects, nCases, nControls, nSNPs)
  names(outputData) <- c("OutSubjects", "OutCases", "OutControls", "OutSNPs")
  # return nSubjects, nMale, nFemale, nSNPs
  return(c(initData, outputData))
}


>>>>>>> reportWriter
