.libPaths("/home/gdc/public/Ref/R")
library(dplyr)
library(magrittr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- args[2]

ancestry <- read.table(paste0(dir, "/ancestry_", name, ".txt"), header=T)
sample <- ancestry[,1]
Q_data <- ancestry[,-1]
index <- apply(Q_data, 1, function(x) which(x==max(x)))
ancestry_names <- colnames(Q_data)
ancestry_decision <- data.frame(sample, index)
colnames(ancestry_decision) <- c("ID", "code_number")
ancestry_decision$ancestry <- ancestry_names[index]
ancestry_decision$prediction_percentage <- apply(Q_data, 1, max) * 100

fam_name <- paste0("study.", name, ".unrelated.fam")
fam_path <- paste0(dir, "/relatedness/", fam_name)
fam_file <- read.table(fam_path, header=FALSE)
colnames(fam_file) <- c("FID", "IID", "MID", "PID", "gender", "phenotype")
fam_file$ID <- paste0(fam_file$FID, "_", fam_file$IID)

joined_file <- dplyr::inner_join(fam_file, ancestry_decision, by="ID") %>% 
  dplyr::select(all_of(c("FID", "IID", "ancestry", "prediction_percentage", "gender", "phenotype")))
joined_file$ancestry[which(joined_file$prediction_percentage < 80)] <- "Other"


output_name <- paste0("study.",name,".unrelated.comm.popu")
output_path <- paste0(dir, "/", output_name)
write.table(joined_file, file=output_path, row.names = F, col.names = F, quote=F, sep="\t")


