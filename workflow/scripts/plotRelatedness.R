library(tidyverse)
theme_set(theme_minimal())
args <- commandArgs(trailingOnly = TRUE)
kin_file <- args[1]
kinid_file <- args[2]
out <- args[3]

ids <- read_table(kinid_file, col_names = c("FID", "IID")) |>
  mutate(idx = row_number()) |>
  unite(ID, FID, IID, sep = "_")

lines <- readLines(kin_file)
lower <- lapply(strsplit(lines, "\\s+"), as.numeric)
n <- length(lower)
K <- matrix(NA, n, n)
for (i in seq_len(n)) {
  K[i, 1:i] <- lower[[i]]
}
K[upper.tri(K)] <- t(K)[upper.tri(K)]

kin <- as.data.frame(K) |>
  mutate(ID1 = ids$ID) |>
  pivot_longer(-ID1, names_to = "j", values_to = "KINSHIP") |>
  mutate(j = as.integer(str_remove(j, "V")), ID2 = ids$ID[j]) |>
  filter(ID1 != ID2) |>
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
