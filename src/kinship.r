args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- args[2]

kin <- read.table(paste0(dir,"/",name), header = TRUE)
library(reshape2)
kin_matrix <- dcast(kin, ID1 ~ ID2, value.var = "Kinship")
rownames(kin_matrix) <- kin_matrix\$ID1
kin_matrix <- as.matrix(kin_matrix[, -1])
saveRDS(kin_matrix, file = paste0(dir,"/kinship_matrix.rds")