## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(Hmisc)
library(dendextend)
library(adegenet)
library(fpc)

## loading the gds of the data and pullling some attributes out
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset.gds")
sample.id <- read.gdsn(index.gdsn(wheat.subset, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat.subset, "snp.id"))
snp.pos <- read.gdsn(index.gdsn(wheat.subset, "snp.position"))
snp.chrom <- read.gdsn(index.gdsn(wheat.subset, "snp.chromosome"))
genotypes <- read.gdsn(index.gdsn(wheat.subset, "genotype"))
# meta.data <- cbind(BP = read.gdsn(index.gdsn(wheat.subset, "samp.annot/BP")),
#                    Year = read.gdsn(index.gdsn(wheat.subset, "samp.annot/Year")),
#                    MC = read.gdsn(index.gdsn(wheat.subset, "samp.annot/MC")),
#                    Desig = read.gdsn(index.gdsn(wheat.subset, "samp.annot/designation")),
#                    Origin = read.gdsn(index.gdsn(wheat.subset, "samp.annot/origin")),
#                    Strength = read.gdsn(index.gdsn(wheat.subset, "samp.annot/strength")),
#                    Colour = read.gdsn(index.gdsn(wheat.subset, "samp.annot/colour")),
#                    Season = read.gdsn(index.gdsn(wheat.subset, "samp.annot/season")))
snpgdsClose(wheat.subset)
# ## Revalue some data
# meta.data$Year <- cut(as.numeric(meta.data$Year), 
#                       breaks = c(1800, 1900, 1920, 1940, 1960, 1980, 2000, 2020), labels = F)
# meta.data$Year[is.na(meta.data$Year)] <- 8
# revaluing and imputing gentotypes
genotypes <- replace(genotypes, genotypes == 3, NA)
genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes2 <- knncatimpute(genotypes)

## finding the cluster order
#creating the dist object
# wheat.dist <- as.dist(1-cor(genotypes2, use = "pairwise.complete.obs"))
# # computing likeliest pamk groups
# pop.code.pamk <- factor(pamk(wheat.dist, krange=2:10, diss = T)$pamobject$clustering)

## pop.code from dendrogram clustering
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\sample_indices.Rdata")
dend <- genotypes.pvclust %>% as.dendrogram
pop.code.dend <- vector()
pop.code.dend[labels(dend)] <- c(rep("HRS", 238), rep("Chinese", 10), rep("Australian", 4),
                                 rep("PSR", 28), rep("White Spring", 47), 
                                 rep("Winter - Out", 3), rep("Winter", 58))
pop.code.dend <- pop.code.dend[-sample.indices]

## restore genotype categories
genotypes2 <- replace(genotypes2, genotypes2 == 1, 0)

## making the genlight object
wheat.subset <- new("genlight", gen = t(data.matrix(genotypes2)), ploidy = 2, 
                  ind.names = sample.id, loc.names = snp.id, 
                  chromosome = snp.chrom, position = snp.pos,
                  pop = pop.code.dend)

## saving the genlight object
save(wheat.subset, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genlight.RData")