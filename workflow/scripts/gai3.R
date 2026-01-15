library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- args[2]

ancestry <- list.files(path = paste0(dir, "/rfmix"), pattern = "ancestry_chr.*.rfmix.Q", full.names=T) |>
  map(read_table, skip = 1) |>
  reduce(left_join, by = "#sample") |>
  # make one set of columns appear with the intended ancestry names
  rename_with(
    .fn = ~ str_replace(., pattern = "\\.[x]$", replacement = ""),
  ) |>
  # take the mean across all chromosomes
  mutate(
    `#sample`,
    across(
    .cols = !starts_with("#sample") & !matches("\\.x$|\\.y$"),
    .fns = list(mean = ~ mean(c_across(starts_with(str_sub(cur_column(), 1,3))), na.rm = T)),
    .names = "{str_sub(.col, 1, 3)}"
  )) |>
  select(-c(ends_with(".x"), ends_with(".y")))

sample <- ancestry[,1]
Q_data <- ancestry[,-1]
index <- unlist(apply(Q_data, 1, function(x) which.max(x), simplify = T))
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


output_path <- paste0(dir, "/ancestry_estimation/study.",name,".unrelated.comm.popu")
write.table(joined_file, file=output_path, row.names = F, col.names = F, quote=F, sep="\t")


