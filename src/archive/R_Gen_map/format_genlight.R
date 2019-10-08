## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(Hmisc)
library(dendextend)
library(adegenet)
library(fpc)

# adegenetTutorial("genomics")

## indices
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\snp_indices.RData") 

## loading the gds of the data and pullling some attributes out
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data_2.gds")
pop.code.bp <- read.gdsn(index.gdsn(wheat.data, "BP"))
pop.code.or <- read.gdsn(index.gdsn(wheat.data, "origin"))
pop.code.st <- read.gdsn(index.gdsn(wheat.data, "strength"))
pop.code.co <- read.gdsn(index.gdsn(wheat.data, "colour"))
pop.code.se <- read.gdsn(index.gdsn(wheat.data, "season"))
pop.code.de <- read.gdsn(index.gdsn(wheat.data, "designation"))
pop.code.mc <- read.gdsn(index.gdsn(wheat.data, "MC"))
sample.id <- read.gdsn(index.gdsn(wheat.data, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
snp.pos <- read.gdsn(index.gdsn(wheat.data, "snp.position"))
snp.chrom <- read.gdsn(index.gdsn(wheat.data, "snp.chromosome"))
genotypes <- read.gdsn(index.gdsn(wheat.data, "genotype"))
## Binning the year of release pop groups
# pop.code.year <- read.gdsn(index.gdsn(wheat.data, "Year"))
# pop.code.era <- as.numeric(cut2(as.numeric(as.character(pop.code.year)), g=5))
# pop.code.era[which(is.na(pop.code.era))] <- "UNKNOWN"
# pop.code.era <- factor(pop.code.era)
pop.code.era <- cut(as.numeric(as.character(read.gdsn(index.gdsn(wheat.data, "Year")))),
                breaks = c(1800, 1900, 1920, 1940, 1960, 1980, 2000, 2020), labels = F)
pop.code.era[is.na(pop.code.era)] <- 8


## finding the cluster order
pamk <- pamk(dist(1-abs(cor(genotypes[kept.indices,]))), krange=2:15, diss = T)

## pop.code from dendrogram clustering
pop.code.pamk <- factor(pamk$pamobject$clustering)

## pop.code from dendrogram clustering
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
dend <- genotypes.pvclust %>% as.dendrogram
pop.code.dend <- vector()
pop.code.dend[labels(dend)] <- c(rep("HRS", 238), rep("Chinese", 10), rep("Australian", 4),
                                 rep("PSR", 28), rep("White Spring", 47), 
                                 rep("Winter - Out", 3), rep("Winter", 58))

# Removing SNPs with high MR and low MAF
snpset <- snpgdsSelectSNP(wheat.data, autosome.only = F, maf = 0.05, missing.rate = 0.1)
snpset.id <- unlist(snpset)
informative <- match(snpset.id, snp.id)
genotypes <- read.gdsn(index.gdsn(wheat.data, "genotype"))[informative,]
genotypes <- replace(genotypes, genotypes == 3, NA)
snpgdsClose(wheat.data)

## making the genlight object
wheat.data <- new("genlight", gen = t(data.matrix(genotypes)), ploidy = 2, 
            ind.names = sample.id, loc.names = snp.id[informative], 
            chromosome = as.character(snp.chrom)[informative], position = as.character(snp.pos)[informative],
            pop = pop.code.dend)

# other = list("bp" = pop.code.bp, "era" = pop.code.era, "origin" = pop.code.or, 
#              "strength" = pop.code.se, "colour" = pop.code.co, "season" = pop.code.se,
#              "designation" = pop.code.de, "mc" = pop.code.mc, "dend" = pop.code.dend)

## saving the genlight object
save(wheat.data, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genlight.RData")