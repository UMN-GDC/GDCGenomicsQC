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
msp_list <- list.files(pattern = ".msp.tsv")
numfiles <- length(Qfiles)

print("Input ancestry codes")
ancestry_code <- system(paste("head -n 1 ", msp_list[1]), intern=T)
ancestry_code <- sub(".*:", "", ancestry_code)
ancestry_code <- gsub("\t", ",", ancestry_code)
ancestry_code <- unlist(strsplit(ancestry_code, ","))
ancestry_code <- gsub("\\s+", "", ancestry_code)
code <- sub("=.*", "", ancestry_code)
number <- sub(".*=", "", ancestry_code)

qnames <- c("ID", code)
temp <- read.table(Qfiles[1])
ancestry_array <- array(0, dim=c(dim(temp)[[1]], dim(temp)[[2]], length(code)))

for(i in 1:length(code)){
  temp <- read.table(Qfiles[1])
  colnames(temp) <- qnames
  ancestry_array[,,i] <- as.matrix(temp)
}

sample_names <- temp$ID
rm(temp)

ancestry_array <- ancestry_array[,-1,]
ancestry_array <- array(as.numeric(ancestry_array), dim = dim(ancestry_array))
ancestry_mat <- apply(ancestry_array, c(1,2), mean)
colnames(ancestry_mat) <- code
rownames(ancestry_mat) <- sample_names

print("Outputs")
output_name <- paste0("ancestry_",name,".txt")
output_path <- paste0(dir, "/", output_name)
write.table(ancestry_mat, file=output_path, row.names = T, col.names = T, quote=F)


index <- apply(ancestry_mat, 1, function(x) which(x==max(x)))
#ancestry_vec <- unlist(lapply(index, function(x) return(names(x))))
ancestry_vec <- unname(index)
id_vec <- names(index)
ancestry_decision <- data.frame(id_vec, ancestry_vec)
colnames(ancestry_decision) <- c("ID", "code_number")
ancestry_decision$ancestry <- ancestry_decision$code_number
for(i in 1:length(code)){
  ancestry_decision$ancestry[which(ancestry_decision$ancestry==i)] <- code[i]
}
# Add prediction percentage
ancestry_decision$prediction_percentage <- apply(ancestry_mat, 1, max) * 100

fam_name <- paste0("study.", name, ".unrelated.fam")
fam_path <- paste0(dir, "/relatedness/", fam_name)
fam_file <- read.table(fam_path, header=FALSE)
colnames(fam_file) <- c("FID", "IID", "MID", "PID", "gender", "phenotype")
fam_file$ID <- paste0(fam_file$FID, "_", fam_file$IID)

joined_file <- dplyr::inner_join(fam_file, ancestry_decision, by="ID") %>% 
  dplyr::select(all_of(c("FID", "IID", "ancestry", "prediction_percentage", "gender", "phenotype")))


output_name <- paste0("study.",name,".unrelated.comm.popu")
output_path <- paste0(dir, "/", output_name)
write.table(joined_file, file=output_path, row.names = F, col.names = F, quote=F, sep="\t")



# ancestry_asignment <- matrix(0, nrow=dim(ancestry_mat)[[1]], ncol=dim(ancestry_mat)[[2]])
# colnames(ancestry_asignment) <- code
# rownames(ancestry_asignment) <- sample_names
