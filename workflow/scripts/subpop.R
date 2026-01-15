args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
name <- args[2]

work_dir <- paste0(dir, "/ancestry_estimation")
setwd(work_dir)

x <- read.table("subpop.txt")
data <- read.table("data.txt")

for (i in 1:nrow(x)) {
  y <- data[data[,3]== x[i,2], c(1,2)]
  write.table(y, file = paste(x[i,2]), sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
