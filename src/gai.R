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

# Qfiles <- list.files(pattern = ".rfmix.Q")
# temp <- read.table("/projects/standard/gdc/and02709/ancestry/rff/ancestry_chr1.rfmix.Q")
# qnames <- system(paste("head -n 2 ", Qfiles[1]), intern=T)
# qnames <- qnames[2]
# qnames <- sub("#", "", qnames)
# qnames <- gsub("\t", ",", qnames)
# qnames <- unlist(strsplit(qnames, ","))
# colnames(temp) <- qnames

# Load the RFMix output file
print("Load RFMix output files")
msp_list <- list.files(pattern = ".msp.tsv")
numfiles <- length(msp_list)

print("Initialize null data frame")
msp_data <- data.frame()

print("Input ancestry codes")
ancestry_code <- system(paste("head -n 1 ", msp_list[1]), intern=T)
ancestry_code <- sub(".*:", "", ancestry_code)
ancestry_code <- gsub("\t", ",", ancestry_code)
ancestry_code <- unlist(strsplit(ancestry_code, ","))
ancestry_code <- gsub("\\s+", "", ancestry_code)
code <- sub("=.*", "", ancestry_code)
number <- sub(".*=", "", ancestry_code)

print("Naming msp columns")
msp_colnames <- system(paste("head -n 2 ", msp_list[1]), intern=T)
msp_colnames <- msp_colnames[2]
msp_colnames <- sub("#", "", msp_colnames)
msp_colnames <- gsub("\t", ",", msp_colnames)
msp_colnames <- unlist(strsplit(msp_colnames, ","))

print("Table of outputs")
for(i in 1:numfiles){
  temp <- read.table(msp_list[1], header = F, sep = "\t")
  msp_data <- rbind(msp_data, temp)
}
#msp_data <- read.table(msp_file, header = TRUE, sep = "\t")
rm(temp)
print("Rename msp data")
colnames(msp_data) <- msp_colnames

# Calculate the segment length
print("Extract segment lengths")
msp_data <- msp_data %>%
  mutate(sl = epos - spos)

segment_data <- msp_data %>% dplyr::select(c("chm", "spos", "epos", "sl", "sgpos", "egpos", "n snps"))
colnames(segment_data) <- c("chm", "spos", "epos", "slength", "sgpos", "egpos", "n_snps")
msp_data <- msp_data %>% dplyr::select(-c("chm", "spos", "epos", "sl", "sgpos", "egpos", "n snps"))

print("Initialize ancestry array")
ancestry_array <- array(0, dim=c(dim(msp_data)[[1]], dim(msp_data)[[2]], length(code)), 
                        dimnames=list(rownames(msp_data), colnames(msp_data), code))

print("Fill in ancestry array")
for(i in 1:length(code)){
  index <- which(msp_data==number[i], arr.ind=T)
  index <- cbind(index, rep(i, dim(index)[[1]]))
  ancestry_array[index] <- 1
}

print("Calculations")
length_vec <- segment_data$slength
length_mat <- t(as.matrix(length_vec))
length_tot <- sum(length_vec)
prop_mat <- length_mat/length_tot
hap_mat <- matrix(0, nrow=dim(ancestry_array)[[2]], ncol=dim(ancestry_array)[[3]])
rownames(hap_mat) <- colnames(msp_data)
colnames(hap_mat) <- code
for(i in 1:length(code)){
  hap_mat[,i] <- prop_mat %*% ancestry_array[,,i]
}

nsamples <- dim(hap_mat)[[1]]/2

#ancestry_mat <- matrix(0, nrow=nsamples, ncol=length(code))
#sample_names <- c()

print("haplotype results")
result_list <- lapply(1:nsamples, function(x){
  (hap_mat[2*x-1,]+hap_mat[2*x,])/2
})

sample_names <- str_extract(rownames(hap_mat), "^[^.]+")[2*(1:nsamples)]
#names(result_list) <- sample_names

ancestry_mat <- t(data.frame(result_list))

rownames(ancestry_mat) <- sample_names
colnames(ancestry_mat) <- code

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


