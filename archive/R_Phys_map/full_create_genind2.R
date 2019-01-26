## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(adegenet)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_subset_sample.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
snpgdsClose(wheat)

# making the data palatable by genind
genotypes <- replace(genotypes, genotypes == 0, "A")
genotypes <- replace(genotypes, genotypes == 2, "B")
genotypes <- replace(genotypes, genotypes == 3, "N")

# get the pop code
# legend=c("Hard Red Spring 1", "Hard Red Spring 2", "Soft White Spring", 
#          "Prairie Spring Red", "Hard Red Spring 3", "Hard Red Winter")
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\full_kmeans.RDATA")
pop.code <- bests

## for making subsets of the samples for the different categories excluding unknown samples
index.hrs <- c(which(pop.code == 1), which(pop.code == 2), which(pop.code == 5))
index.sws <- which(pop.code == 3)
index.hrw <- which(pop.code == 6)

## strata
strata.hrs_sws <- data.frame(c(rep("HRS", length(index.hrs)), rep("SWS", length(index.sws))))
colnames(strata.hrs_sws) <- "hrs_sws"
strata.hrs_hrw <- data.frame(c(rep("HRS", length(index.hrs)), rep("HRW", length(index.hrw))))
colnames(strata.hrs_hrw) <- "hrs_hrw"


## making the genind objects
wheat.poppr.hrs_sws <- df2genind(t(data.frame(genotypes[,c(index.hrs, index.sws)])), 
                                 ind.names = as.character(sample.id)[c(index.hrs, index.sws)], 
                                 loc.names = snp.id, NA.char = "N", ploidy = 1, type = "codom", 
                                 ncode = 1, strata = strata.hrs_sws)
wheat.poppr.hrs_hrw <- df2genind(t(data.frame(genotypes[,c(index.hrs, index.hrw)])), 
                                 ind.names = as.character(sample.id)[c(index.hrs, index.hrw)], 
                                 loc.names = snp.id, NA.char = "N", ploidy = 1, type = "codom", 
                                 ncode = 1, strata = strata.hrs_hrw)

# saving the genlight object
save(wheat.poppr.hrs_sws, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\genind_hrs_sws2.RData")
save(wheat.poppr.hrs_hrw, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\genind_hrs_hrw2.RData")
