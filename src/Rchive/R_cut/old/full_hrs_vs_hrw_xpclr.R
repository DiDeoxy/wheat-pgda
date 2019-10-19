library(SNPRelate)
# install.packages("gdata")
library(gdata)
library(plyr)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")


wheat <-
  snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")

## making .geno files
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
desig <- read.gdsn(index.gdsn(wheat, "samp.annot/designation"))

index.hrs <- which(desig == "HRS")
index.hrw <- which(desig == "HRW")
index.hrw <- which(desig == "SWS")

genotypes[genotypes == 3] <- 9
genotypes <- t(interleave(genotypes, genotypes))

write.table(genotypes[which(desig == "HRS"),], file = "Data\\Intermediate\\XPCLR\\HRS.geno", quote = F, row.names = F, col.names = F)
write.table(genotypes[which(desig == "HRW"),], file = "Data\\Intermediate\\XPCLR\\HRW.geno", quote = F, row.names = F, col.names = F)
write.table(genotypes[which(desig == "SWS"),], file = "Data\\Intermediate\\XPCLR\\SWS.geno", quote = F, row.names = F, col.names = F)

## making .snp file
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
snp.chr <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))

data <- read.csv("Data\\Raw\\Genotypes\\Jan_6_wheat_genotypes_curtis.csv",
                 header = T, comment.char = "", quote="", stringsAsFactors = F, row.names = 2)
colnames(data)[1:6] <- data[2,1:6]
data <- data[-2,]

snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snpgdsClose(wheat)

alleles <- list()
# ?data.frame
for (i in 1:13192) {
  first <- sample(c("A", "C", "G", "T"), 1)
  if (first == "A") {
    second <- sample(c("C", "G", "T"), 1)
  } else if (first == "C") {
    second <- sample(c("A", "G", "T"), 1)
  } else if (first == "G") {
    second <- sample(c("A", "C", "T"), 1)
  } else {
    second <- sample(c("C", "G", "T"), 1)
  }
  alleles[[i]] <- cbind(first, second)
}
alleles <- do.call(rbind, alleles)

snps <- cbind(snp.id, snp.chr, data$Position[match(snp.id, row.names(data))], snp.pos, alleles)

write.table(snps, "Data\\Intermediate\\XPCLR\\wheat.snps", quote = F, row.names = F, col.names = F)
