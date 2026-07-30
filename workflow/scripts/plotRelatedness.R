library(tidyverse)
theme_set(theme_minimal())
args <- commandArgs(trailingOnly = TRUE)
kin_file <- args[1]
out <- args[2]

lines <- readLines(kin_file)
kinship <- unlist(lapply(seq_along(lines), function(i) {
  vals <- as.numeric(strsplit(lines[i], "\\s+")[[1]])
  vals[1:(i - 1)]
}))

expected_lines <- data.frame(
  kinship = c(0.25, 0.125, 0.0625, 0.03125),
  label = c("1st-degree", "2nd-degree", "3rd-degree", "4th-degree")
)

tibble(KINSHIP = kinship) |>
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
