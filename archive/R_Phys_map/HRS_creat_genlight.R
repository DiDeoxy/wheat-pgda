## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(dendextend)
library(adegenet)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_imputed_subset.gds")
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chrom <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
genotypes <- replace(genotypes, genotypes == 3, NA)
snpgdsClose(wheat)

# neighbour joining
neighbour <- nj(dist.gene(genotypes, method = "percentage"))
labels(neighbour) <- sample.id
neighbour <- root(neighbour, which(sample.id == "Selkirk"), r = TRUE)
neighbour <- chronos(neighbour, lambda = 0, model = "correlated")
neighbour <- as.dendrogram(neighbour)

pop.code <- vector()
pop.code[match(labels(neighbour), sample.id)] <- c(rep("Hard Red Spring", 204), rep("Prairie Spring Red", 23),
                                                   rep("Hard Red Spring 2", 9), rep("Out", 1))
pop.code <- factor(pop.code)

## making the genlight object
wheat.genlight <- new("genlight", gen = data.matrix(genotypes), ploidy = 2, 
                      ind.names = sample.id, loc.names = snp.id, 
                      chromosome = snp.chrom, position = snp.pos,
                      pop = pop.code)

## saving the genlight object
save(wheat.genlight, file = 
       "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_imputed_subset_genlight.RData")