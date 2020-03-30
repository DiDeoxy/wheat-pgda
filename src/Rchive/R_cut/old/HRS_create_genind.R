## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(adegenet)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("Data\\Formatted\\HRS_phys_subset_sample.gds")
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snpgdsClose(wheat)

# making the data palatable by genind
genotypes <- replace(genotypes, genotypes == 0, "A")
genotypes <- replace(genotypes, genotypes == 2, "B")
genotypes <- replace(genotypes, genotypes == 3, "N")

## Getting the pop code
load(file = "Data\\Formatted\\hrs_kmeans.RDATA")
## getting the kmeans goups
index.hrs2_hrs3 <- c(which(bests == 2), which(bests == 4))
index.hrs1_psr <- c(which(bests == 1), which(bests == 3))

## strata
# strata.psr_hrs <- replace(bests, bests == 2, 1)
# strata.psr_hrs <- data.frame(replace(strata.psr_hrs, strata.psr_hrs == 4, 1))
strata.hrs1_psr <- data.frame(bests[index.hrs1_psr])
strata.hrs2_hrs3 <- data.frame(bests[index.hrs2_hrs3])
# colnames(strata.psr_hrs) <- "psr_hrs"
colnames(strata.hrs1_psr) <- "psr_hrs"
colnames(strata.hrs2_hrs3) <- "hrs2_hrs3"

## making the genind objects
# wheat.poppr.psr_hrs <- df2genind(t(data.frame(genotypes[,index.psr_hrs])), 
#                              ind.names = as.character(sample.id)[index.psr_hrs], 
#                              loc.names = snp.id, NA.char = "N", ploidy = 1,
#                              type = "codom", ncode = 1, strata = strata.psr_hrs)
wheat.poppr.psr_hrs <- df2genind(t(data.frame(genotypes[,index.hrs1_psr])), 
                                 ind.names = as.character(sample.id)[index.hrs1_psr], 
                                 loc.names = snp.id, NA.char = "N", ploidy = 1,
                                 type = "codom", ncode = 1, strata = strata.hrs1_psr)
wheat.poppr.hrs2_hrs3 <- df2genind(t(data.frame(genotypes[,index.hrs2_hrs3])),
                             ind.names = as.character(sample.id)[index.hrs2_hrs3], 
                             loc.names = snp.id, NA.char = "N", ploidy = 1,
                             type = "codom", ncode = 1, strata = strata.hrs2_hrs3)

# saving the genlight object
save(wheat.poppr.psr_hrs, file = "Data\\Formatted\\genind_psr_hrs.RData")
save(wheat.poppr.hrs2_hrs3, file = "Data\\Formatted\\genind_hrs2_hrs3.RData")
