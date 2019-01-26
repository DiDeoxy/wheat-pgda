## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(adegenet)
library(poppr)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H.DESKTOP-AJ57KB6\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))

# making the data palatable by genind
genotypes <- replace(genotypes, genotypes == 0, "A")
genotypes <- replace(genotypes, genotypes == 2, "B")
genotypes <- replace(genotypes, genotypes == 3, "N")


pop.code.era <- cut(as.numeric(as.character(read.gdsn(index.gdsn(wheat, "samp.annot/Year")))),
                    breaks = c(1800, 1900, 1920, 1940, 1960, 1980, 2000, 2020), labels = F)
snpgdsClose(wheat)

## for making subsets of the samples for the different categories excluding unknown samples
# index.1 <- which(pop.code.era == 1)
# index.2 <- which(pop.code.era == 2)
# index.3 <- which(pop.code.era == 3)
index.4 <- which(pop.code.era == 4)
index.5 <- which(pop.code.era == 5)
index.6 <- which(pop.code.era == 6)
index.7 <- which(pop.code.era == 7)

## strata
strata.45 <- as.data.frame(pop.code.era[c(index.4, index.5)])
colnames(strata.45) <- "four_five"
strata.56 <- as.data.frame(pop.code.era[c(index.5, index.6)])
colnames(strata.56) <- "five_six"
strata.67 <- as.data.frame(pop.code.era[c(index.6, index.7)])
colnames(strata.67) <- "six_seven"
strata.47 <- as.data.frame(pop.code.era[c(index.4, index.7)])
colnames(strata.47) <- "four_seven"

## making the genind objects
wheat.poppr.45 <- df2genind(t(data.frame(genotypes[,c(index.4, index.5)])), 
                            ind.names = sample.id[c(index.4, index.5)], 
                            loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.45)
wheat.poppr.56 <- df2genind(t(data.frame(genotypes[,c(index.5, index.6)])), 
                            ind.names = sample.id[c(index.5, index.6)], 
                            loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.56)
wheat.poppr.67 <- df2genind(t(data.frame(genotypes[,c(index.6, index.7)])), 
                            ind.names = sample.id[c(index.6, index.7)], 
                            loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.67)
wheat.poppr.47 <- df2genind(t(data.frame(genotypes[,c(index.4, index.7)])), 
                            ind.names = sample.id[c(index.4, index.7)], 
                            loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.47)

# amova
wheat.amova <- poppr.amova(wheat.poppr.45, hier = ~four_five, missing = "genotype", within = F, clonecorrect = F)
wheat.amova

wheat.amova <- poppr.amova(wheat.poppr.56, hier = ~five_six, missing = "genotype", within = F, clonecorrect = F)
wheat.amova

wheat.amova <- poppr.amova(wheat.poppr.67, hier = ~six_seven, missing = "genotype", within = F, clonecorrect = F)
wheat.amova

wheat.amova <- poppr.amova(wheat.poppr.47, hier = ~four_seven, missing = "genotype", within = F, clonecorrect = F)
wheat.amova
