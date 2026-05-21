library(tidyverse)
theme_set(theme_minimal())
args <- commandArgs(trailingOnly = TRUE)
hwe_file <- args[1]
out <- args[2]

hwe <- read_table(hwe_file, comment = "#")
p_col <- if ("P" %in% colnames(hwe)) "P" else tail(colnames(hwe), 1)

hwe |>
  mutate(log10p = -log10(.data[[p_col]])) |>
  filter(is.finite(log10p)) |>
  ggplot(aes(x = log10p)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = -log10(1e-6), color = "red", linetype = "dashed") +
  geom_vline(xintercept = -log10(1e-10), color = "darkred", linetype = "dashed") +
  labs(x = expression(-log[10](HWE ~ p)), y = "Count",
       caption = "Red dashed lines at p = 1e-6 and p = 1e-10")
ggsave(out, width = 9, height = 5)
