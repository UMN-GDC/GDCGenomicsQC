args <- commandArgs(trailingOnly = TRUE)
het_file <- args[1]
output_dir <- args[2]

het <- read.table(het_file, head=FALSE, col.names=c("FID", "IID", "O_HOM", "E_HOM", "OBS_CT", "F"))
het_fail = subset(het, (het$F < mean(het$F)-3*sd(het$F)) | (het$F > mean(het$F)+3*sd(het$F)));
if (nrow(het_fail) > 0) {
    het_fail$F_DST = (het_fail$F - mean(het$F)) / sd(het$F);
    write.table(het_fail[, c("FID", "IID")], paste0(output_dir, "/het_fail_ind.txt"), row.names=FALSE)
}
