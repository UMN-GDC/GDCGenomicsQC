#' @title Read the initial data from the log file'
#' @description Read the initial data from the log file considering expanded plink options'
#' @param filepath The path to the log file'
#' @return A list of named vectors characterizing the input and output data'
#' @export
#' @examples
#' extractLog("../QCtutorial/logs/wgas2.log")

suppressMessages(library(tidyverse)) 

extractLog_v2<- function(filepath, plinkoption){
  file <- file(paste0(wd,filepath), "r")
  log <- readLines(file)
  close(file)
  
  
  str_1 <- log[grep("people .* from .fam", log)]
  test1 <- str_extract_all(str_1, "[0-9]+")
  nSubjects <- as.numeric(test1[[1]][1])
  nMale <- as.numeric(test1[[1]][2])
  nFemale <- as.numeric(test1[[1]][3])
  
  # SNPsNSubjects()
  str_2 <- log[grep("loaded from .bim", log)]
  nSNPs <- as.numeric(gsub("[^0-9]", "", str_2))
  
  # Phenotype values
  str_3 <- log[grep("phenotype values loaded from .fam", log)]
  nPheno <- as.numeric(gsub("[^0-9]", "", str_3))
  
  initData <-cbind(nSubjects, nMale, nFemale, nSNPs, nPheno)
  colnames(initData) <- c("InSubjects", "InMale", "InFemale", "InSNPs", "InPheno")
  
  # 5x2 table
  
  # Passing QC
  # SNPs and subejcts
  
  # #Initializing indicators
  # out_label="NotUsed"
  # num_args=1
  
  # Section for implementing switch to look for different things  
  if(length(plinkoption) == 0) {plinkoption=2} 
  if (plinkoption == 1) { #--indep-pairwise
    string_to_find_a="Pruned .* variants from chromosome .* leaving" 
    # out_label="NumPeopleRemoved"
  }
  if (plinkoption==2) { #Default
    error_message= "Invalid plink option selected"
    stop(print(error_message))
  }
  
  string_to_find_b="[0-9]+ variants and [0-9]+ people pass filters and QC." #All logs have this in common!
  # For plink option 1
  string_to_find_c="Pruning complete.* variants removed"
  
  # Extracting specific to plink option provided
  ## Untested
  str_out_unique<- str_extract_all(log[grep(string_to_find_a, log)], "[0-9]+")
  if(length(str_out_unique) == 0) {
    error_message= "Invalid log file provided for this plink option selected"
    stop(print(error_message))
  } else{
    
    num_rows=22 # Hard coded for testing purposes
    df = data.frame(matrix( 
      vector(), num_rows, 3, dimnames=list(c(), c("PrunedSNPS","Chr","RemainingSNPS"))), 
      stringsAsFactors=F) 
    
    for(i in 1:num_rows) {
      current_var <- paste0("row_", i)
      current_row <-as.numeric(str_out_unique[[i]])
      df$PrunedSNPS[i] <- current_row[1]
      df$Chr[i] <- current_row[2]
      df$RemainingSNPS[i] <- current_row[3]
    }
    
    # row1 <- as.numeric(str_out_unique[[1]])
    # test_Chr <- as.numeric(str_out_unique[[2]])
    # test_Remaining <- as.numeric(str_out_unique[[3]])
    # test_table=cbind(test_Pruned,test_Chr,test_Remaining)
    # colnames(test_table) <- c("PrunedSNPS", "Chr", "RemainingSNPS")
    # 3x23 table
  }
  
  #### Lots of small alterations need to be done below with the change from 'c' to 'cbind' in the function ####
  
  # Retreiving surviving #SNPs and participants
  str_out_1<- str_extract_all(log[grep(string_to_find_b, log)], "[0-9]+")
  if(length(str_out_1)==0){
    error_message_2="Bad string_to_find_b"
    stop(print(error_message_2))
  }

  
  # Part 2 of log output
  str_out_2 <- str_extract_all(log[grep(string_to_find_c, log)], "[0-9]+")
  if(length(str_out_2)==0){
    error_message_2="Bad string_to_find_c"
    stop(print(error_message_2))
  }
  nSNPSCut <- as.numeric(str_out_2[[1]][1])
  nSNPSTot <- as.numeric(str_out_2[[1]][2])

  # Part 1 of log output
  nSNPS2 <- as.numeric(str_out_1[[1]][1])
  nSubjs2 <- as.numeric(str_out_1[[1]][2])
  outputData <- cbind(nSubjs2, nSNPSCut, nSNPS2)
  
  # Allows for different column names between plink options
  # outputData <- cbind(....)
  colnames(outputData) <- c("OutSubjects", "NumSNPStoPrune", "OutSNPs")
  # 3x2 table
  
  # Part 3 of log output
  
  # Merging information
  nrow1=nrow(initData)
  nrow2=nrow(outputData)
  if(nrow1 > nrow2){
    error_message_3="Number of rows for initData are greater than the outputData"
    stop(print(error_message_3, nrow1, nrow2))
  }
  if(nrow2 > nrow1){
    error_message_3="Number of rows for outputData are greater than the initData"
    stop(print(error_message_3, nrow2, nrow1))
  }

  
  final_output_1 = cbind(initData, outputData)
  
  # return nSubjects, nMale, nFemale, nSNPs
  return(list(x=final_output_1, y=df))
}


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#### NOT fixed /adjusted yet ####
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
final_name_2 =paste0(wd, "/each_SNP_", output_name)
print(final_name)

plink_selected <- switch(plink_option,
                         "indep-pairwise"=1,
                         "hwe" = 6)

print(plink_selected)
table1=extractLog_v2(filename, plink_selected)
table1$x #Shows the number of SNPS removed?
print(names(table1))
table2 <- table1$y
table1$y #Useful table 
print(" ")
write.table(table1$x, file = paste0(final_name), quote = F, row.names = F, col.names = T)
write.table(table2, file = paste0(final_name_2), quote = F, row.names = F, col.names = T)
#After extractLog need to save it as a table to be called within quarto


