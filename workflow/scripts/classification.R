library(argparse)
library(tidyverse) |> suppressPackageStartupMessages()
library(randomForest) |> suppressPackageStartupMessages()

parser <- ArgumentParser(description= "Estimate ancestry off of multiple predictors.")
parser$add_argument("--eigen_ref", type= "character", default = NULL,
    help = "Filepath to the eigenvector file (plink style)")
parser$add_argument("--eigen_sample", type= "character", default = NULL,
    help = "Filepath to the eigenvector file (plink style)")
parser$add_argument("--npc", type= "integer", default = NULL,
    help = "Number of PCs to use for nacestry prediction. Does not apply to VAE and UMAP predictions.")
parser$add_argument("--labels", type= "character", default = NULL,
    help = "Path to labels file. Assumes file columns, FID, IID, Population")
parser$add_argument("--vae", type= "character", default=NULL,
    help = "Filepath to the vae file (popvae style)")
parser$add_argument("--umap_ref", type= "character", default= NULL,
    help = "Filepath to the umap file")
parser$add_argument("--umap_sample", type= "character", default= NULL,
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

ref <- read_table(args$labels) |>
  select(FID = FamilyID, IID = SampleID, POP = Superpopulation)

print("this is the unique")
print(unique(ref$POP))

set.seed(args$seed)

# load dimension reductions
PCs <- read_table(args$eigen_ref, col_names= TRUE) |> select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM))
colnames(PCs) <- c("IID", paste0("pc_", 1:(ncol(PCs) -1) ))
ref <- full_join(ref, PCs, by = c("IID")) |> drop_na(pc_1)
pcMod <- randomForest::randomForest(formula = factor(POP) ~ pc_1 + pc_2 + pc_3 + pc_4 + pc_5 + pc_6+ pc_7 + pc_8 + pc_9 + pc_10, data = ref) 
saveRDS(pcMod, paste0(args$out, "/RFpc.Rds"))
ref <- ref |> mutate(pc_label = pcMod$predicted)

sampleDF <- read_table(args$eigen_sample, col_names= TRUE) |> select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM))
colnames(sampleDF) <- c("IID", paste0("pc_", 1:(ncol(sampleDF) -1) ))
sampleDF$pc_label <- predict(pcMod, sampleDF)

if (!is.null(args$umap_ref)) {
  umap_ref <- read_csv(args$umap_ref)
  colnames(umap_ref) <- c("IID", str_replace(colnames(umap_ref)[-c(1)], "UMAP", "umap_"))
  ref <- full_join(ref, umap_ref, by = c("IID")) |> drop_na(umap_1)
  umapMod <- randomForest::randomForest(formula = factor(POP) ~ umap_1 + umap_2, data = ref)
  saveRDS(umapMod, paste0(args$out, "/RFumap.Rds"))
  ref$umap_label <- predict(umapMod, ref)

  umap_sample <- read_csv(args$umap_sample)
  colnames(umap_sample) <- c("IID", str_replace(colnames(umap_sample)[-c(1)], "UMAP", "umap_"))
  sampleDF <- full_join(sampleDF, umap_sample, by = c("IID")) |> drop_na(umap_1)
  sampleDF$umap_label <- predict(umapMod, sampleDF)

}
if (!is.null(args$vae)) {
  vae <- read_table(args$vae, col_names= TRUE)
  colnames(vae) <- paste0("vae_", colnames(vae))
  colnames(vae)[length(colnames(vae))] <- "IID"
  ref <- full_join(ref, vae, by = c("IID")) |> drop_na()
  vaeMod <- randomForest::randomForest(formula = factor(POP) ~ vae_mean1 + vae_mean2, data = ref)
  saveRDS(vaeMod, paste0(args$out, "/RFvae.Rds"))
  ref <- ref |> mutate(vae_label = vaeMod$predicted)

  # Need to do VAE on ref and then project sample into those VAEs
}

backgroundDat <- ref |>
   select(-any_of(c("vae_sd1", "vae_sd2", "PAT", "MAT", "PHENO"))) |>
   pivot_longer(
     -c(IID, POP),
     #names_pattern = "([A-Za-z]+?)(\\d+|label)$",
     names_pattern = "([a-z]{2,4})_([a-z]*[0-9]{0,2})$",
     names_to = c("alg", ".value"),
   ) |>
  # mutate(
  #   `1` = coalesce(`1`, mean1),
  #   `2` = coalesce(`2`, mean2),
  # ) |>
  select(-c(starts_with("mean")))
foregroundDat <- sampleDF |>
   select(IID, starts_with("pc_"), starts_with("umap_")) |>
   pivot_longer(
     -c(IID),
     #names_pattern = "([A-Za-z]+?)(\\d+|label)$",
     names_pattern = "([a-z]{2,4})_([a-z]*[0-9]{0,2})$",
     names_to = c("alg", ".value"),
   ) |>
  # mutate(
  #   `1` = coalesce(`1`, mean1),
  #   `2` = coalesce(`2`, mean2),
  # ) |>
  select(-c(starts_with("mean")))


backgroundDat |>
  filter(alg=="umap") |>
  ggplot(aes(x = `1`, y = `2`, fill = label, group = interaction(label, alg))) + 
  stat_density_2d(geom = "polygon", alpha =0.25, aes(fill = label)) +
  theme_minimal() +
  geom_point(data = foregroundDat |> filter(alg=="umap"),
  aes(x = `1`, y = `2`, fill = label, ), shape = 21 , col = "black") +
  ylab("UMAP embedding 2") +
  xlab("UMAP embedding 1")
ggsave(paste0(args$out, "/UMAP_referenceSpace.svg"), dpi = 300, width = 9, height= 5)
backgroundDat |>
  filter(alg=="pc") |>
  ggplot(aes(x = `1`, y = `2`, fill = label, col = label, group = interaction(label, alg))) + 
  #stat_density_2d(geom = "polygon", alpha =0.25, aes(fill = label)) +
  geom_point(alpha =0.25) +
  theme_minimal() +
  geom_point(data = foregroundDat |>
    filter(alg=="pc"),
  aes(x = `1`, y = `2`, fill = label, ), shape = 21 , col = "black") +
  ylab("PC 2") +
  xlab("PC 1")
ggsave(paste0(args$out, "/PC_referenceSpace.svg"), dpi = 300, width = 9, height= 5)

ref |>
  ggplot(aes(x = pc_1, pc_2)) +
  geom_point()

sampleDF |>
  relocate(IID) |>
  write_delim(paste0(args$out, "/latentDistantRelatedness.tsv"), delim = "\t")
