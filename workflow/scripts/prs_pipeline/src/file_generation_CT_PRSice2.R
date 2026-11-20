#!/usr/bin/env Rscript
#Adding in ability for script to accept command line arguments

## Inputs
# fam.file:    study sample fam file
# qassoc:      training_sumstats_corrected.assoc.linear
# bim.file:    SNP information (.bim file) from population that summary statistics were generated from
# out.path:    Path for outputs to be stored

args = commandArgs(trailingOnly=TRUE)

# test if all arguments are provided: if not, return an error
if (length(args)!= 4) {
  stop("All four arguments must be provided. 'qassoc.path bim.file.path table.with.beta.se.path out.path'.", call.=FALSE)
}

args[1] -> qassoc.path # Full path to the a summary stat file that has these columns TEST, SNP, & BP
args[2] -> bim.file.path # SNP information (.bim file) from population that summary statistics were generated from
args[3] -> table.with.beta.se.path # Summary stat file with a header that contains columns rsid, chr, & beta_se
args[4] -> out.path # Path for outputs to be stored

## load necessary libraries
suppressMessages(library(dplyr))

# generate summary stats file so that it is in the correct format
## read in summary statistics file (same code that was used in LDpred2/lassosum2 code implementation)
qassoc <- read.table(qassoc.path, header=TRUE)
qassoc <- qassoc %>% filter(TEST == "ADD")
qassoc <- qassoc %>% rename(rsid = SNP, POS = BP)
bim.file <- bigreadr::fread2(bim.file.path, select = c(1, 4, 5, 6)) # reads in relevant information only (not the entire file; check to see if this needs to be the GWAS or study sample here??)
colnames(bim.file) <- c("chr", "pos", "bim.a1", "bim.a0")
table.with.beta.se <- read.table(table.with.beta.se.path, header = TRUE, sep = "", stringsAsFactors = FALSE)
beta.se <- table.with.beta.se %>%  select(rsid, chr, beta_se)
## merge bim, beta.se, and qassoc together
qassoc <- qassoc %>% rename(chr = CHR, pos = POS)
sumstats <- qassoc %>% transmute(rsid, chr, pos, a1=A1, beta=BETA, n_eff=NMISS) %>% inner_join(bim.file, by=c("chr", "pos")) 
sumstats <- sumstats %>% mutate(a0=ifelse(a1==bim.a1, bim.a0, bim.a1))
sumstats <- sumstats %>% inner_join(beta.se, by=c("chr", "rsid"))
## rename columns so it matches the expected format
sumstats.renamed <- sumstats %>% rename(CHR = chr, SNP= rsid,A1 = a1, A2 = a0, BP = pos)
## remove duplicated columns so that there are no issues
sumstats.renamed <- sumstats.renamed %>% select(SNP, CHR, BP, A1, A2, beta, beta_se, n_eff)
## rearrange columns so that it matches the expected format
sumstats.rearranged <- sumstats.renamed %>% relocate(SNP, CHR, BP, A1, A2, beta, beta_se, n_eff) 
## add p-value column
z <- sumstats.rearranged$beta / sumstats.rearranged$beta_se # make z scores
p <- 2*(1-pnorm(abs(z)))
sumstats.rearranged$P <- p
sumstats.rearranged <- sumstats.rearranged %>% relocate(P, .after = beta_se) 

summary.stats.output=file.path(out.path, "CT_PRSice2_summary_stat_file.txt")
write.table(sumstats.rearranged, 
            file =summary.stats.output, 
            sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
head(sumstats.rearranged)