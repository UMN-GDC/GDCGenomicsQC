library(tidyverse)
theme_set(theme_minimal())
args <- commandArgs(trailingOnly = TRUE)

smissFiles <- unlist(strsplit(args[1], " "))
vmissFiles <- unlist(strsplit(args[2], " "))
imissOut <- args[3]
lmissOut <- args[4]

smiss <- lapply(smissFiles, read_table) |> bind_rows()
vmiss <- lapply(vmissFiles, read_table) |> bind_rows()

smiss |>
  ggplot(aes(x = F_MISS * 100)) +
  geom_histogram() +
  geom_vline(xintercept = c(0.1, 0.02), color = "red") +
  xlab("Subject genome missingness (%)")
ggsave(imissOut, width = 9, height= 5)
  
vmiss |>
  ggplot(aes(x = F_MISS * 100)) +
  geom_histogram() +
  geom_vline(xintercept = c(0.1, 0.02), color = "red") +
  xlab("SNP missingness (%)")
ggsave(lmissOut, width = 9, height= 5)
