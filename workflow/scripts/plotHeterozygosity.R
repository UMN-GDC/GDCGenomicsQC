library(tidyverse)
theme_set(theme_minimal())
args <- commandArgs(trailingOnly = TRUE)
het_file <- args[1]
out <- args[2]

het <- read_table(het_file)
het |>
  ggplot(aes(x = F)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean(het$F) + c(-3, 3) * sd(het$F),
             color = "red", linetype = "dashed") +
  labs(x = "Inbreeding coefficient (F)", y = "Count",
       caption = "Red dashed lines at mean ± 3 SD")
ggsave(out, width = 9, height = 5)
