.libPaths("/projects/standard/gdc/public/Ref/R")
library(dplyr)
library(magrittr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- args[2]

print("Initialize rfmix files")
work_dir <- paste0(dir, "/rfmix")
setwd(work_dir)
Qfiles <- list.files(pattern = ".rfmix.Q")
temp <- read.table(paste0(work_dir,"/ancestry_chr1.rfmix.Q"))
qnames <- system(paste("head -n 2 ", Qfiles[1]), intern=T)
qnames <- qnames[2]
qnames <- sub("#", "", qnames)
qnames <- gsub("\t", ",", qnames)
qnames <- unlist(strsplit(qnames, ","))
colnames(temp) <- qnames
sample <- temp[,1]

# Load the RFMix output file
print("Load RFMix Q files")
Q_list <- lapply(Qfiles, function(f) {
  dat <- read.table(f, header=FALSE)
  colnames(dat) <- qnames
  return(dat)
})

print("Converting Q_list into 3D array")
# Dimensions
num_chromosomes <- length(Q_list)
num_individuals <- nrow(Q_list[[1]])
num_ancestries <- ncol(Q_list[[1]])

# Initialize array
Q_array <- array(0, dim = c(num_individuals, num_ancestries, num_chromosomes),
                 dimnames = list(NULL, qnames, paste0("chr", 1:num_chromosomes)))

# Fill array with numeric ancestry values
for (i in 1:num_chromosomes) {
  Q_array[,,i] <- apply(Q_list[[i]], 2, as.numeric)
}

Q_data <- apply(Q_array, c(1, 2), mean)
ancestry_mat <- data.frame(Q_data)
ancestry_mat$sample <- sample


print("Outputs")
output_name <- paste0("ancestry_",name,".txt")
output_path <- paste0(dir, "/", output_name)
write.table(ancestry_mat, file=output_path, row.names = T, col.names = T, quote=F)
Q_data <- Q_data[,-1]

index <- apply(Q_data, 1, function(x) which(x==max(x)))
#ancestry_vec <- unlist(lapply(index, function(x) return(names(x))))
ancestry_names <- qnames[-1]
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


