# Creates a single report from the QC steps

#These two files are created and overwrote whenever there is a --missing call
# For Longs QC this is at the end of QC steps

#Adding in ability for script to accept command line arguments
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least two arguments: if not, return an error
if (length(args)<=1) {
  stop("A path to the output location and the file need to be provided.", call.=FALSE)
} else if (length(args)>=2) {
  args[1] -> base_path #Place where the tables are located
  path_to_logs=paste0(base_path, "/data/")
  # default output file
  args[2] -> file_preffix 
  if(length(args)==2){
    output_location <- "results"
    path_to_save_report=paste0(base_path, "/results/")
  }
  if(length(args)>2) {
    args[3] -> output_location
    path_to_save_report=paste0(base_path, output_location)
  }

}

wd=path_to_logs
print(wd)
print(file_preffix)
print(output_location)


# read data into R from temporary space
indmiss<-read.table(file="plink.imiss", header=TRUE)
snpmiss<-read.table(file="plink.lmiss", header=TRUE)
gender <- read.table("gender_check.sexcheck", header=T,as.is=T)
maf_freq <- read.table("MAF_check.frq", header =TRUE, as.is=T)
hwe<-read.table (file="plink.hwe", header=TRUE)
hwe_zoom<-read.table (file="plinkzoomhwe.hwe", header=TRUE)
het <- read.table("R_check.het", head=TRUE)
relatedness = read.table("pihat_min0.2.genome", header=T)
relatedness_zoom = read.table("zoom_pihat.genome", header=T)

# print("Checking tables are read properly")
# head(indmiss)
# head(snpmiss)
# head(gender)
# head(maf_freq)
# head(hwe)
# head(hwe_zoom)
# head(het)
# head(relatedness)
# head(relatedness_zoom)


setwd("/panfs/jay/groups/16/saonli/baron063/R")

suppressMessages(library(tidyverse))

#Back to where the log data is
setwd(wd)
print("Files available to be selected in data directory")
list.files()

#Reading in the tables for later use
QC2_geno_table <- read.table(file = "QC2_geno.txt",  quote="\"", comment.char="")

QC3_mind_table <- read.table("QC3_mind.txt",  quote="\"", comment.char="")
QC4_geno_table <- read.table("QC4_geno.txt",  quote="\"", comment.char="")
QC5_mind_table <- read.table("QC5_mind.txt",  quote="\"", comment.char="")
QC6_sex_check_table <- read.table("QC6_sex_check.txt", quote="\"", comment.char="")
QC7_maf_table <- read.table("QC_7_maf.txt", quote="\"", comment.char="")
QC8_hwe_table <- read.table("QC_8_hwe.txt", quote="\"", comment.char="")
QC8b_hwe_table <- read.table("QC_8b_hwe.txt", quote="\"", comment.char="")
QC9_filter_founders_table <- read.table("QC_9_filter-founders.txt", quote="\"", comment.char="")

QC_indep_pairwise_table <- read.csv("QC_indep_pairwise.txt", sep = "")
QC_indep_pairwise_bychr <- read.csv("each_SNP_QC_indep_pairwise.txt", sep = "")

#Key 
# common            InSubjects InMale InFemale InSNPs 
# common            OutSubjects *unique* OutSNPs
#
# Start of unique element based on plink option chosen
# --mind            NumPeopleRemoved
# --geno            NumVariantsRemoved
# --check-sex       Numx NumY NumProblems
# --maf             NumVariantsRemoved
# --filter-founders NumRemoved NumFounder NumNonFounder
# --hwe             NumVariantsRemoved
#
# Below is from the logReader_extended function
# --indep-pairwise  NumSNPStoPrune <-- table1
#                   PrunedSNPS Chr RemainingSNPS <-- table2

#Going to where the report should be saved
setwd(path_to_save_report)

#### Start of QCreport_2.pdf ####
pdf("QCreport_2.pdf")
geno_table_summary = rbind(QC2_geno_table, QC4_geno_table)
# Probably will turn it into a gt?
# Would be good to rearrange the columns? 
# These all have 7 columns
mind_table_summary = rbind(QC3_mind_table, QC5_mind_table)

## Testing
grid.table(mind_table_summary, show.rownames = FALSE)

QC6_sex_check_table
# Turn it into a gt?

QC7_maf_table
# gt?

hwe_check_table = rbind(QC8_hwe_table, QC8b_hwe_table)
# gt?

QC9_filter_founders_table

QC_indep_pairwise_table

QC_indep_pairwise_bychr
# gt?

dev.off() # Ends pdf creation. 

#### Start of QCreport.pdf ####

pdf("QCreport.pdf") #indicates pdf format and gives title to file

print(length(indmiss))

print(nrow(indmiss))

#Option 1
data.frame("Subject" = 1:nrow(indmiss),
          "Missingness" = indmiss$F_MISS) %>%
  ggplot(aes(x = Missingness)) +
  geom_histogram() + 
  geom_vline(xintercept = 0.15, color = "red") + 
  xlab("Missingness per subject") +
  ggtitle("% SNPS missing per subject")
#Option 2
hist(indmiss[,6],main="Histogram individual missingness", xlab = "Proportion of missing SNPs") #selects column 6, names header of file

#Option 1
data.frame("SNP" = 1:nrow(snpmiss),
          "Missingness" = snpmiss$F_MISS) %>%
  ggplot(aes(x = Missingness)) +
  geom_histogram() +
  geom_vline(xintercept = 0.15, color = "red") +
  xlab("Missingness per SNP") +
  ggtitle("% calls missing per SNP")
#Option 2
hist(snpmiss[,5],main="Histogram SNP missingness", xlab = "Proportion of sample missing for each SNP")  

# print("hist_miss.R Script Success!")

# chromosome homozygosity estimate (F statistic) is < 0.2 for Females and as males if the estimate is > 0.8

#### Base R versions 
## Not yet implemented into the Quarto document...
hist(gender[,6],main="Gender", xlab="F Value")

male=subset(gender, gender$PEDSEX==1)
hist(male[,6],main="Men",xlab="F Value")

female=subset(gender, gender$PEDSEX==2)
hist(female[,6],main="Women",xlab="F Value")

## Adding a barplot to show the number that are problems 
temptab= table(gender$STATUS)
barplot(temptab, main = "Homozygosity Analysis", xlab = "Status")

# print("gender_check.R Script Success!")
#Option 2 MAF
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
#####

#Option 1 MAF
data.frame("SNP" = 1:nrow(maf_freq),
           "MAF" = maf_freq[,5]) %>%
  ggplot(aes(x = MAF)) +
  geom_histogram() +
  geom_vline(xintercept = 0.05, color = "red") +
  xlab("Minor allele frequency") +
  ggtitle("MAF distribution")

# print("MAF_check.R Script Success!")

#Option 2 HWE
hist(hwe[,9],main="Histogram HWE", xlab = "P-value")
#Option 1 HWE 
data.frame("SNP" = 1:nrow(hwe),
          "HWE" = hwe[,9]) %>%
          ggplot(aes(x = HWE)) +
          geom_histogram() +
          geom_vline(xintercept = -0.05, color = "red") +
          xlab("log(HWE p-value)") +
          ggtitle("Hardy-Weinberg Equilibrium p-value distribution")

#Option 2 b HWE
hist(hwe_zoom[,9],main="Histogram HWE: strongly deviating SNPs only", xlab = "P-value")

# print("hwe.R Script Success!")

#Option 1 Heterozygosity
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
data.frame("Subject" = 1:nrow(het),
          "Heterozygosity" = het$HET_RATE) %>%
  ggplot(aes(x = Heterozygosity)) +
  geom_histogram() + 
  geom_vline(xintercept = c(-2, 2), color = "red") + 
  xlab("Heterozygosity F statistic") +
  ggtitle("Heterozygous F statistic distribution")
#Option 2 Heterozygosity
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")

#Option 2 Heterozygosity further analysis 
# Should make a quarto / ggplot version 
# ## Adding a barplot to show the number that are problems 
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));

placeholder = c()
for(i in 1:nrow(het)) {
  if(het[i, 2] %in% het_fail[[2]]) {
    placeholder[i]="PROBLEM"
  }
  if(!het[i,2] %in% het_fail[[2]]) {
    placeholder[i] = "OK"
  }
}

temp_table = table(placeholder)
barplot(temp_table, main = "Heterozygosity Analysis", xlab = "Status")

# print("check_heterozygosity_rate.R Script Success!")

#Option 2 Relatedness plots
## Should I include these?
### If so I'd need to make them ggplots
par(pch=16, cex=1)
plot(relatedness$Z0, relatedness$Z1, xlim=c(0,1), ylim=c(0,1), xlab = "P(IBD=0)", 
     ylab = "P(IBD=1)", main = "Relatedness")

plot(density(relatedness$Z0), xlab ="P(IBD=0)", main = "Relatedness Density Plot P(IBD=0)")

plot(density(relatedness$Z1), xlab ="P(IBD=1)", main = "Relatedness Density Plot P(IBD=1)")

with(subset(relatedness,RT=="PO") , points(Z0,Z1,col=4))
with(subset(relatedness,RT=="UN") , points(Z0,Z1,col=3))

table(relatedness$Z0)
table(relatedness$Z1)

par(pch=16, cex=1)
plot(relatedness_zoom$Z0, relatedness_zoom$Z1, xlim=c(0,0.02), ylim=c(0.98,1),
     xlab = "P(IBD=0)", ylab = "P(IBD=1)", main = "Relatedness Zoom")
with(subset(relatedness_zoom,RT=="PO") , points(Z0,Z1,col=4))
with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col=3))

hist(relatedness[,10],main="Histogram relatedness", xlab= "Proportion IBD")  

# print("Relatedness.R Script Success!")

# print("Unified Report Completed!")


dev.off() # Ends pdf creation. 


