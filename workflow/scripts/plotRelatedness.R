library(tidyverse)
theme_set(theme_minimal())
args <- commandArgs(trailingOnly = TRUE)
kin_file <- args[1]
kinid_file <- args[2]
out <- args[3]

lines <- readLines(kin_file)
lower <- lapply(strsplit(lines, "\\s+"), as.numeric)
n <- length(lower)

id_lines <- readLines(kinid_file)
ids <- read.table(text = id_lines[1:n], col.names = c("FID", "IID")) |>
  unite(ID, FID, IID, sep = "_")

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

expected_lines <- data.frame(
  kinship = c(0.25, 0.125, 0.0625, 0.03125),
  label = c("1st-degree", "2nd-degree", "3rd-degree", "4th-degree")
)

kin |>
  ggplot(aes(x = KINSHIP)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(data = expected_lines, aes(xintercept = kinship),
             color = "gray50", linetype = "dotted") +
  geom_text(data = expected_lines, aes(x = kinship, label = label),
            y = 0, hjust = -0.1, vjust = 1.5, size = 3, color = "gray50") +
  geom_vline(xintercept = 0.0442, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.0884, color = "darkred", linetype = "dashed") +
  geom_vline(xintercept = 0.354, color = "orange", linetype = "dashed") +
  annotate("text", x = 0.05, y = 0, label = "2nd-degree", hjust = 0, vjust = 0, size = 3, color = "red") +
  annotate("text", x = 0.09, y = 0, label = "3rd-degree", hjust = 0, vjust = 0, size = 3, color = "darkred") +
  labs(x = "KING kinship coefficient", y = "Count",
       caption = "Dotted: expected kinship | Dashed: filtering thresholds")
ggsave(out, width = 9, height = 5)
