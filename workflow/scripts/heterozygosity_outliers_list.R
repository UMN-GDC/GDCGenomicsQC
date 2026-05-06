args <- commandArgs(trailingOnly = TRUE)
het_file <- args[1]
output_dir <- args[2]

het <- read.table(het_file, head=FALSE, col.names=c("FID", "IID", "O_HOM", "E_HOM", "N_NM"))
het$HET_RATE = (het$N_NM - het$O_HOM) / het$N_NM
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
if (nrow(het_fail) > 0) {
    het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
    write.table(het_fail[, c("FID", "IID")], paste0(output_dir, "/het_fail_ind.txt"), row.names=FALSE)
}
