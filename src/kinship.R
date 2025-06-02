args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- args[2]

kin <- read.table(paste0(dir,"/",name), header = TRUE)
library(reshape2)
kin_matrix <- dcast(kin, ID1 ~ ID2, value.var = "Kinship")
rownames(kin_matrix) <- kin_matrix$ID1
kin_matrix <- as.matrix(kin_matrix[, -1])
saveRDS(kin_matrix, file = paste0(dir,"/kinship_matrix.rds"))
coords <- which(kin_matrix > 0.05, arr.ind = TRUE)
row_ids <- rownames(kin_matrix)[coords[,1]]
col_ids <- colnames(kin_matrix)[coords[,1]]
ids_all <- c(row_ids, col_ids)
ids_unique <- unique(ids_all)
df <- data.frame(ids_unique, ids_unique)
write.table(df, file="to_exclude.txt", col.names = F, row.names = F, quote=F)