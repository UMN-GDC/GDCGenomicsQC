#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(jsonlite)
})

parser <- ArgumentParser(description = "Create GRM-ordered visit phenotypes for AdjHE")
parser$add_argument("--longitudinal", required = TRUE)
parser$add_argument("--grm-id", required = TRUE)
parser$add_argument("--out", required = TRUE)
parser$add_argument("--manifest", required = TRUE)
args <- parser$parse_args()

long <- fread(args$longitudinal)
grm <- fread(args$grm_id, header = FALSE)
setnames(grm, c("FID", "IID"))
grm[, `:=`(FID = as.character(FID), IID = as.character(IID), grm_order = .I)]
long[, `:=`(FID = as.character(FID), IID = as.character(IID))]

visits <- sort(unique(long$visit))
wide <- dcast(long, FID + IID ~ visit, value.var = "pheno")
visit_cols <- setdiff(names(wide), c("FID", "IID"))
setnames(wide, visit_cols, paste0("visit_", visits))
visit_cols <- paste0("visit_", visits)

wide <- merge(grm, wide, by = c("FID", "IID"), all.x = TRUE, sort = FALSE)
setorder(wide, grm_order)
wide[, grm_order := NULL]
wide <- wide[complete.cases(wide[, ..visit_cols])]
if (!nrow(wide)) stop("No complete visit phenotypes remained after GRM-ID matching")

fwrite(wide, args$out, sep = " ", col.names = FALSE)
write_json(
  list(
    visits = as.list(visits),
    phenotype_names = as.list(visit_cols),
    n_final = nrow(wide)
  ),
  args$manifest,
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("Wrote visit phenotype file:", args$out, "\n")
cat("Subjects:", nrow(wide), "visits:", length(visit_cols), "\n")
