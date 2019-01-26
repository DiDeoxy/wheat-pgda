## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(dendextend)
library(adegenet)
# library(fpc)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_subset_both.gds")
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chrom <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
genotypes <- replace(genotypes, genotypes == 3, NA)
snpgdsClose(wheat)

pop.code <- bests

## making the genlight object
wheat.genlight <- new("genlight", gen = data.matrix(t(genotypes)), ploidy = 2, 
                    ind.names = sample.id, loc.names = snp.id, 
                    chromosome = snp.chrom, position = snp.pos,
                    pop = pop.code)

## saving the genlight object
save(wheat.genlight, file = 
       "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_subset_both_genlight.RData")
