library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()
library(uwot)

parser <- ArgumentParser(description= "Develop a UMAP embedding off of PCA compressed genomics data.")
parser$add_argument("--eigens", type= "character", 
    help = "Filepath to the eigenvector file (plink style)")
parser$add_argument("--sample", type= "character", 
    help = "Filepath to the subject scores")
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
samplePCs <- data.table::fread(args$sample)
#samplePCs <- data.table::fread("../../toyPipeline/01-globalAncestry/sampleRefPCscores.sscore") |> 

if (is.null(args$npc)) {
  npc <- ncol(pcs) -4
} else if (args$npc > (ncol(pcs) -4)) {
  print("Desired number of PCs exceeds dimension of eigenvec file. Selecting all PCs.")
  npc <- ncol(pcs) - 4
} else {
  npc <- args$npc
}

mod <- pcs |>
  select(starts_with("PC")) |>
  scale() |>
  as.data.frame() |>
  umap(
    n_threads = args$threads,
    n_components = args$ncoords,
    n_neighbors = args$neighbors, ret_mod = T)
    # n_threads = 1,
    # n_components = 2,
    # n_neighbors = 50, ret_mod = T)
studyUmap <- pcs |>
  select(starts_with("PC")) |>
  scale() |>
  as.data.frame() |>
  umap_transform(model = mod)


mod$embedding |>
  magrittr::set_colnames(paste("UMAP", 1:args$ncoords, sep = "")) |>
  cbind(pcs[,c("#FID", "IID")]) |>
  relocate(`#FID`, IID) |>
  data.table::fwrite(file = paste0(args$out, "_ref.csv"),
         row.names = FALSE)
print(paste0("UMAP coordinates saved to ", args$out))

studyUmap |>
  magrittr::set_colnames(paste("UMAP", 1:args$ncoords, sep = "")) |>
  cbind(samplePCs[, c("#FID", "IID")]) |>
  relocate(`#FID`, IID) |>
  data.table::fwrite(file = paste0(args$out, "_sample.csv"),
         row.names = FALSE)
