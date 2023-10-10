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


