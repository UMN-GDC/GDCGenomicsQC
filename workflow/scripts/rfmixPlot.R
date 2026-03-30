library(GenomicRanges)
library(tidyverse)
library(vroom)
library(RIdeogram)

genomeBuild <- "hg38"
mspPath <- "../../toyPipeline/02-localAncestry/chr20.lai.msp.tsv"

data(human_karyotype, package= "RIdeogram")

subj <- "SAS34825_SAS34825"
hap0 <- paste0(subj, ".0")
hap1 <- paste0(subj, ".1")
# header <- names(read_tsv(mspPath, n_max=0, show_col_types=FALSE, skip=1))
colors <- c(
  "0" = "#E69F00",
  "1" = "#56B4E9",
  "2" = "#009E73",
  "3" = "#009E73"
)

msp <- vroom(mspPath, show_col_types = FALSE,
  col_names = TRUE, skip =1,
  col_sel = c(`#chm`, spos, epos, sgpos, egpos,
  !!hap0, !!hap1)) |>
  rename(Chr = `#chm`, Start = spos, End = epos, Value = !!hap0, hap2 = !!hap1)

p = ideogram(karyotype = human_karyotype, overlaid = msp, label_type = "marker", colorset1 = colors )
convertSVG(p, "chr20", device="png")

