## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(plyr)
library(Hmisc)
library(dendextend)
library(adegenet)
library(fpc)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_imputed_sample_subset.gds")
pop.code.st <- read.gdsn(index.gdsn(wheat, "samp.annot/strength"))
pop.code.st <- replace(pop.code.st, pop.code.st == "UNKNOWN", "N/A")
pop.code.co <- read.gdsn(index.gdsn(wheat, "samp.annot/colour"))
pop.code.co <- replace(pop.code.co, pop.code.co == "UNKNOWN", "N/A")
pop.code.se <- read.gdsn(index.gdsn(wheat, "samp.annot/season"))
pop.code.se <- replace(pop.code.se, pop.code.se == "UNKNOWN", "N/A")
snpgdsClose(wheat)

## for making subsets of the samples for the different categories excluding unknown samples
index.st
index.co
index.se
index.hr_sw <- c(which(pop.code.st == "Hard")[which(pop.code.st == "Hard") %in% which(pop.code.co == "Red")],
                 which(pop.code.st == "Soft")[which(pop.code.st == "Soft") %in% which(pop.code.co == "White")])
c(index.hr, index.sw)

## strata
strata.st <- data.frame(pop.code.st[which(pop.code.st != "N/A")])
colnames(strata.st) <- "strength"
strata.co <- data.frame(pop.code.co[which(pop.code.co != "N/A")])
colnames(strata.co) <- "colour"
strata.se <- data.frame(pop.code.se[which(pop.code.se != "N/A")])
colnames(strata.se) <- "season"
strata.hr_sw <- data.frame(pop.code.st)
colnames(strata.hr_sw) <- "hr_sw"

## the actual data
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))

# making the data palatable by genind
genotypes <- replace(genotypes, genotypes == 0, "A")
genotypes <- replace(genotypes, genotypes == 2, "B")
genotypes <- replace(genotypes, genotypes == 3, "N")

## making the genind objects
wheat.poppr.st <- df2genind(t(data.frame(genotypes[,index.st])), 
                            ind.names = as.character(sample.id)[index.st], 
                            loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                            ncode = 1, strata = strata.st)
wheat.poppr.co <- df2genind(t(data.frame(genotypes[,index.co])), 
                            ind.names = as.character(sample.id)[index.co], 
                            loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                            ncode = 1, strata = strata.co)
wheat.poppr.se <- df2genind(t(data.frame(genotypes[,index.se])), 
                            ind.names = as.character(sample.id)[index.se], 
                            loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                            ncode = 1, strata = strata.se)
wheat.poppr.hr_sw <- df2genind(t(data.frame(genotypes)), 
                               ind.names = as.character(sample.id), 
                               loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                               ncode = 1, strata = strata.hr_sw)

# saving the genlight object
save(wheat.poppr.st, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_st.RData")
save(wheat.poppr.co, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_co.RData")
save(wheat.poppr.se, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_se.RData")
save(wheat.poppr.hr_sw, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_hr_sw.RData")
