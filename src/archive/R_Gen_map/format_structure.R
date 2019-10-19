## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(Hmisc)
library(dendextend)

load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\snp_indices.RData") ## indices

## loading the gds of the data and pullling some attributes out
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data.gds")
genotypes <- read.gdsn(index.gdsn(wheat.data, "genotype"))
pop.code.mc <- read.gdsn(index.gdsn(wheat.data, "MC"))
pop.code.bp <- read.gdsn(index.gdsn(wheat.data, "BP"))
sample.id <- read.gdsn(index.gdsn(wheat.data, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
snp.pos <- read.gdsn(index.gdsn(wheat.data, "snp.position"))

## Binning the year of release pop groups
pop.code.year <- read.gdsn(index.gdsn(wheat.data, "Year"))
pop.code.era <- as.numeric(cut2(as.numeric(as.character(pop.code.year)), g=5))
pop.code.era[which(is.na(pop.code.era))] <- "UNKNOWN"
pop.code.era <- factor(pop.code.era)

## pop.code from dendrogram clustering
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
dend <- as.dendrogram(genotypes.pvclust)
pop.code.dend <- vector()
pop.code.dend[labels(dend)] <- c(rep("RW", 58), rep("WS", 62), rep("RS", 268))
pop.code.dend <- as.factor(pop.code.dend)

## reading in metadata
metadata <- read.csv("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Cultivar Info\\Parsed\\metadata.csv",
                     header = T, stringsAsFactors = T)

## final format
markers <- rbind(snp.id[kept.indices], snp.pos[kept.indices]/1000000)
data <- cbind(as.character(metadata$Genotype.Name), as.numeric(pop.code.mc),
              t(data.matrix(genotypes[kept.indices,])))
# as.numeric(pop.code.bp), as.numeric(pop.code.era), as.numeric(pop.code.dend)

## Writing out
write.table(markers, "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes.txt",
            col.names = F, row.names = F, sep = " ", quote = F)
write.table(data, "C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes.txt",
            col.names = F, row.names = F, sep = " ", quote = F, append = T)