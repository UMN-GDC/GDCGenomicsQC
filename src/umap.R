library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description= "Develop a UMAP embedding off of PCA compressed genomics data.")
parser$add_argument("--eigens", type= "character", 
    help = "Filepath to the eigenvector file (plink style)")
parser$add_argument("--out", type= "character", 
    help = "Filepath to save the projection onto UMAP.")
parser$add_argument("--npc", type= "integer", default = NULL,
    help = "Number of PCs for genotype compression (default is max number from eigenvec file)")
parser$add_argument("--neighbors", type= "integer", default = 50,
    help = "Number of nearest neighbors to define local distance manifold. (default = 50)")
parser$add_argument("--threads", type= "integer", default = 1,
    help = "Number of computation threads. (default = 1)")
parser$add_argument("--ncoords", type= "integer", default = 3,
    help = "Number of embedding coordinates. (default = 3)") # try 2-10
parser$add_argument("--seed", type= "integer", default = as.integer(Sys.time()), 
    help = "Specify the desired heritability. Default system time")
args <- parser$parse_args()

set.seed(args$seed)

pcs <- data.table::fread(args$eigens)

if (is.null(args$npc)) {
  npc <- ncol(pcs) -2
} else if (args$npc > (ncol(pcs) -2)) {
  print("Desired number of PCs exceeds dimension of eigenvec file. Selecting all PCs.")
  npc <- ncol(pcs) -2
} else {
  npc <- args$npc
}

pcs |>
  magrittr::set_colnames(c("FID", "IID", paste("PC", 1:npc, sep = ""))) |>
  scale() |>
  uwot::umap(n_threads = args$threads,
             n_components = args$ncoords,
             n_neighbors = args$neighbors) |>
  #with(layout) |>
  magrittr::set_colnames(paste("UMAP", 1:args$ncoords, sep = "")) |>
  data.table::fwrite(file = paste0(args$out, ".csv"),
         row.names = TRUE)
print(paste0("UMAP coordinates saved to ", args$out, ".csv"))
