# Load necessary library
suppressMessages(library(data.table))

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("A path to the .sscore file needs to be provided.", call.=FALSE)
} else if (length(args)>1) {
    args[1] -> sscore_file_path #Full path to the .sscore file
    args[2] -> output_location #Place to store R2 results
} else {
  args[1] -> sscore_file_path #Full path to the .sscore file
  getwd() -> output_location #Place to store R2 results
}

# Read in the .sscore file
prs_data <- fread(sscore_file_path)

# View the first few rows
head(prs_data)

# Simple model: phenotype ~ PRS
model <- lm(PHENO1 ~ SCORE1_AVG, data = prs_data)

# View R²
result_1=summary(model)$r.squared
final_name_1=file.path(output_location, "PRS_sscore_R_sqr.txt")
write.table(result_1, 
            file =final_name_1, 
            sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)

result_adj=summary(model)$adj.r.squared
final_name_adj=file.path(output_location, "adj_PRS_sscore_Rsqr.txt")
write.table(result_adj, 
            file =final_name_adj, 
            sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
