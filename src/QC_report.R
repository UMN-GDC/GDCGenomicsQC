# Creates a single report from the QC steps

#These two files are created and overwrote whenever there is a --missing call
# For Longs QC this is at the end of QC steps
indmiss<-read.table(file="../data/reportData/plink.imiss", header=TRUE)
snpmiss<-read.table(file="plink.lmiss", header=TRUE)
# read data into R 

pdf("QCreport.pdf") #indicates pdf format and gives title to file
hist(indmiss[,6],main="Histogram individual missingness", xlab = "Proportion of missing SNPs") #selects column 6, names header of file

hist(snpmiss[,5],main="Histogram SNP missingness", xlab = "Proportion of sample missing for each SNP")  

# print("hist_miss.R Script Success!")

# chromosome homozygosity estimate (F statistic) is < 0.2 for Females and as males if the estimate is > 0.8
gender <- read.table("plink.sexcheck", header=T,as.is=T)

hist(gender[,6],main="Gender", xlab="F Value")

male=subset(gender, gender$PEDSEX==1)
hist(male[,6],main="Men",xlab="F Value")

female=subset(gender, gender$PEDSEX==2)
hist(female[,6],main="Women",xlab="F Value")

## Adding a barplot to show the number that are problems 
temptab= table(gender$STATUS)
barplot(temptab, main = "Homozygosity Analysis", xlab = "Status")

# print("gender_check.R Script Success!")

maf_freq <- read.table("MAF_check.frq", header =TRUE, as.is=T)
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")

# print("MAF_check.R Script Success!")


hwe<-read.table (file="plink.hwe", header=TRUE)
hist(hwe[,9],main="Histogram HWE", xlab = "P-value")

hwe_zoom<-read.table (file="plinkzoomhwe.hwe", header=TRUE)
hist(hwe_zoom[,9],main="Histogram HWE: strongly deviating SNPs only", xlab = "P-value")

# print("hwe.R Script Success!")

het <- read.table("R_check.het", head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")


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


relatedness = read.table("pihat_min0.2.genome", header=T)
par(pch=16, cex=1)
plot(relatedness$Z0, relatedness$Z1, xlim=c(0,1), ylim=c(0,1), xlab = "P(IBD=0)", 
     ylab = "P(IBD=1)", main = "Relatedness")

plot(density(relatedness$Z0), xlab ="P(IBD=0)", main = "Relatedness Density Plot P(IBD=0)")

plot(density(relatedness$Z1), xlab ="P(IBD=1)", main = "Relatedness Density Plot P(IBD=1)")

with(subset(relatedness,RT=="PO") , points(Z0,Z1,col=4))
with(subset(relatedness,RT=="UN") , points(Z0,Z1,col=3))

table(relatedness$Z0)
table(relatedness$Z1)

relatedness_zoom = read.table("zoom_pihat.genome", header=T)
par(pch=16, cex=1)
plot(relatedness_zoom$Z0, relatedness_zoom$Z1, xlim=c(0,0.02), ylim=c(0.98,1),
     xlab = "P(IBD=0)", ylab = "P(IBD=1)", main = "Relatedness Zoom")
with(subset(relatedness_zoom,RT=="PO") , points(Z0,Z1,col=4))
with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col=3))

relatedness = read.table("pihat_min0.2.genome", header=T)
hist(relatedness[,10],main="Histogram relatedness", xlab= "Proportion IBD")  

# print("Relatedness.R Script Success!")

# print("Unified Report Completed!")


dev.off() # Ends pdf creation. 


