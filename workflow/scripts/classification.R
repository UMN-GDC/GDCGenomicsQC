library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()
library(randomForest) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description= "Estimate ancestry off of multiple predictors.")
parser$add_argument("--pc", type= "character", default = NULL,
    help = "Filepath to the eigenvector file (plink style)")
parser$add_argument("--npc", type= "integer", default = NULL,
    help = "Number of PCs to use for nacestry prediction. Does not apply to VAE and UMAP predictions.")
parser$add_argument("--labels", type= "character", default = NULL,
    help = "Path to labels file. Assumes file columns, FID, IID, Population")
parser$add_argument("--vae", type= "character", default=NULL,
    help = "Filepath to the vae file (popvae style)")
parser$add_argument("--umap", type= "character", default= NULL,
    help = "Filepath to the umap file")
parser$add_argument("--out", type= "character", default=NULL,
    help = "Filepath to save the projection onto UMAP.")
parser$add_argument("--rseed",
  type= "integer", default = as.integer(Sys.time()),
  help = "Specify the desired heritability. Default system time")
args <- parser$parse_args()
# for development
#args <- parser$parse_args(c("--pc", "testData/1kg.eigenvec",
#  "--umap", "testData/1kgUmap.csv",
#  "--labels", "testData/1000G.GRCh38.popu",
#  "--vae", "testData/1kgvae_latent_coords.txt",
#  "--out", "testData/1kg"
# ))
# dat <- c("testData/1kgThinAFR.fam", "testData/1kgThinEUR.fam") |>
#   map(read_table, col_names = c("FID", "IID", "PAT", "MAT", "POP", "PHENO")) |>
#   reduce(rbind) |>
#   mutate(POP = factor(rep(c("AFR", "EUR"), each = 5000)))

dat <- read_table(args$labels, col_names = c("IID", "POP"))

set.seed(args$seed)

# load dimension reductions
if (!is.null(args$pc)) {
  PCs <- read_table(args$pc, col_names= TRUE)
  colnames(PCs) <- c("FID", "IID", paste0("pc_", 1:(ncol(PCs) -2) ))
  dat <- full_join(dat, PCs, by = c("IID")) |> drop_na()
  pcMod <- randomForest::randomForest(formula = factor(POP) ~ pc_1 + pc_2 + pc_3 + pc_4 + pc_5 + pc_6+ pc_7 + pc_8 + pc_9 + pc_10, data = dat) 
  saveRDS(pcMod, paste0(args$out, "/RFpc.Rds"))
  dat <- dat |> mutate(pc_label = pcMod$predicted)
}
if (!is.null(args$vae)) {
  vae <- read_table(args$vae, col_names= TRUE)
  colnames(vae) <- paste0("vae_", colnames(vae))
  colnames(vae)[length(colnames(vae))] <- "IID"
  dat <- full_join(dat, vae, by = c("IID")) |> drop_na()
  vaeMod <- randomForest::randomForest(formula = factor(POP) ~ vae_mean1 + vae_mean2, data = dat)
  saveRDS(vaeMod, paste0(args$out, "/RFvae.Rds"))
  dat <- dat |> mutate(vae_label = vaeMod$predicted)
}
if (!is.null(args$umap)) {
  umap <- read_csv(args$umap) |>
    select(FID = `#FID`, IID, contains("UMAP"))
  colnames(umap) <- c("FID", "IID", str_replace(colnames(umap)[-c(1,2)], "UMAP", "umap_"))
  dat <- full_join(dat, umap, by = c("FID", "IID")) |> drop_na()
  umapMod <- randomForest::randomForest(formula = factor(POP) ~ umap_1 + umap_2, data = dat)
  saveRDS(umapMod, paste0(args$out, "/RFumap.Rds"))
  dat <- dat |> mutate(umap_label = umapMod$predicted)
}

vizDat <- dat |>
  select(-any_of(c("vae_sd1", "vae_sd2", "PAT", "MAT", "PHENO"))) |>
  pivot_longer(
    -c(FID, IID, POP),
    #names_pattern = "([A-Za-z]+?)(\\d+|label)$",
    names_pattern = "([a-z]{2,4})_([a-z]*[0-9]{0,2})$",
    names_to = c("alg", ".value"),
  ) #|>
  #mutate(
  #  `1` = coalesce(`1`, mean1),
  #  `2` = coalesce(`2`, mean2),
  #) |>
  #select(-c(starts_with("mean")))

vizDat |>
  mutate(Concordant = POP== label) |>
  # mutate(
  #   across(`1`:`20`, scale),
  #   .by = c(alg, mislabeled)
  # ) |>
  ggplot(aes(x= `1`, y = `2`, col = POP, shape = Concordant)) +
  geom_point(alpha = 0.25) +
  facet_wrap(~ alg, scales= "free", labeller = label_both, ncol = 2) +
  # facet_grid(rows = vars(alg), cols = vars(mislabeled), scales= "free", labeller = label_both) +
  theme_minimal() 
ggsave(paste0(args$out, "/latentDistantRelatedness.png"))

dat |>
  relocate(FID, IID) |>
  write_delim(paste0(args$out, "/latentDistantRelatedness.tsv"), delim = "\t")
