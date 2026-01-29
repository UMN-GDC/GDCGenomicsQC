library(tidyverse)
theme_set(theme_minimal())
args <- commandArgs(trailingOnly = TRUE)

imissIn <- args[1]
lmissIn <- args[2]
imissOut <- args[3]
lmissOut <- args[4]

read_table(imissIn) |>
  ggplot(aes(x = F_MISS * 100)) +
  geom_histogram() +
  geom_vline(xintercept = c(0.1, 0.02), color = "red") +
  xlab("Subject genome missingness (%)")
ggsave(imissOut)
  
read_table(lmissIn) |>
  ggplot(aes(x = F_MISS * 100)) +
  geom_histogram() +
  geom_vline(xintercept = c(0.1, 0.02), color = "red") +
  xlab("SNP missingness (%)")
ggsave(lmissOut)
