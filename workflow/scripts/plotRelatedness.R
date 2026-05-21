library(tidyverse)
theme_set(theme_minimal())
args <- commandArgs(trailingOnly = TRUE)
kin_file <- args[1]
out <- args[2]

kin <- read_table(kin_file)
kin <- kin |>
  mutate(pair = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")) |>
  distinct(pair, .keep_all = TRUE)

kin |>
  ggplot(aes(x = KINSHIP)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.0442, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.0884, color = "darkred", linetype = "dashed") +
  geom_vline(xintercept = 0.354, color = "orange", linetype = "dashed") +
  annotate("text", x = 0.05, y = 0, label = "2nd-degree", hjust = 0, vjust = 0, size = 3, color = "red") +
  annotate("text", x = 0.09, y = 0, label = "3rd-degree", hjust = 0, vjust = 0, size = 3, color = "darkred") +
  labs(x = "KING kinship coefficient", y = "Count",
       caption = "Dashed lines at relatedness thresholds: 2nd-degree, 3rd-degree")
ggsave(out, width = 9, height = 5)
