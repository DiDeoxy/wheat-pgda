library(SNPRelate)
library(scrime)
library(BHC)
library(RColorBrewer)
library(extrafont)
library(dendextend)

## loading the data
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset.gds")
genotypes <- t(read.gdsn(index.gdsn(wheat.subset, "genotype")))
sample.id <- as.character(read.gdsn(index.gdsn(wheat.subset, "sample.id")))
snp.pos <- read.gdsn(index.gdsn(wheat.subset, "snp.position"))
snp.chr <- read.gdsn(index.gdsn(wheat.subset, "snp.chromosome"))
meta.data <- cbind(BP = as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/BP"))),
                   Year = as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/Year"))),
                   MC = as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/MC"))),
                   Desig = as.character(read.gdsn(index.gdsn(wheat.subset, "samp.annot/designation"))))
snpgdsClose(wheat.subset)
rownames(meta.data) <- sample.id
## reformat genotypes and impute missing ones
genotypes <- replace(genotypes, genotypes == 3, NA)
genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes2 <- knncatimpute(genotypes)

## bayesian clustering
bhc.dend <- bhc(genotypes2, itemLabels = sample.id, dataType = "multinomial", verbose = T)
save(bhc.dend, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\bhc_dend.Rdata")
