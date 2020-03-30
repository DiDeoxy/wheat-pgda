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
pop.code.bp <- replace(pop.code.bp, pop.code.bp == "UNKNOWN", "N/A")
pop.code.or <- read.gdsn(index.gdsn(wheat.data, "origin"))
pop.code.or <- replace(pop.code.or, pop.code.or == "UNKNOWN", "N/A")
pop.code.st <- read.gdsn(index.gdsn(wheat.data, "strength"))
pop.code.st <- replace(pop.code.st, pop.code.st == "UNKNOWN", "N/A")
pop.code.co <- read.gdsn(index.gdsn(wheat.data, "colour"))
pop.code.co <- replace(pop.code.co, pop.code.co == "UNKNOWN", "N/A")
pop.code.se <- read.gdsn(index.gdsn(wheat.data, "season"))
pop.code.se <- replace(pop.code.se, pop.code.se == "UNKNOWN", "N/A")
pop.code.de <- read.gdsn(index.gdsn(wheat.data, "designation"))
pop.code.de <- replace(pop.code.de, pop.code.de == "UNKNOWN", "N/A")
pop.code.mc <- read.gdsn(index.gdsn(wheat.data, "MC"))
pop.code.mc <- replace(pop.code.mc, pop.code.mc == "UNKNOWN", "N/A")
# Binning the year of release pop groups
pop.code.era <- cut(as.numeric(as.character(read.gdsn(index.gdsn(wheat.data, "Year")))),
    breaks = c(1800, 1900, 1920, 1940, 1960, 1980, 2000, 2020), labels = F)
pop.code.era[which(is.na(pop.code.era))] <- "UNKNOWN"
pop.code.era <- factor(pop.code.era)
# finding the pamk cluster memberships
genotypes <- read.gdsn(index.gdsn(wheat.data, "genotype"))
pamk <- pamk(dist(1-abs(cor(genotypes[kept.indices,]))), krange=2:15, diss = T)
pop.code.pamk <- factor(pamk$pamobject$clustering)

## for making subsets of the samples for the different categories excluding unknown samples
index.bp <- which(pop.code.bp != "N/A")
index.or <- which(pop.code.or != "N/A")
index.st <- which(pop.code.st != "N/A")
index.co <- which(pop.code.co != "N/A")
index.se <- which(pop.code.se != "N/A")
index.de <- which(pop.code.de != "N/A")
index.mc <- which(pop.code.mc != "N/A")
index.hr <- which(pop.code.st == "Hard")[which(pop.code.st == "Hard") %in% which(pop.code.co == "Red")]
index.sw <- which(pop.code.st == "Soft")[which(pop.code.st == "Soft") %in% which(pop.code.co == "White")]
index.hr.sw <- c(index.hr, index.sw)

## strata
strata.bp <- data.frame(pop.code.bp[index.bp])
colnames(strata.bp) <- "bp"
strata.or <- data.frame(pop.code.or[index.or])
colnames(strata.or) <- "origin"
strata.st <- data.frame(pop.code.st[index.st])
colnames(strata.st) <- "strength"
strata.co <- data.frame(pop.code.co[index.co])
colnames(strata.co) <- "colour"
strata.se <- data.frame(pop.code.se[index.se])
colnames(strata.se) <- "season"
strata.de <- data.frame(pop.code.de[index.de])
colnames(strata.de) <- "designation"
strata.mc <- data.frame(pop.code.mc[index.mc])
colnames(strata.mc) <- "mc"
strata.pk <- data.frame(pop.code.pamk)
colnames(strata.pk) <- "pamk"
strata.hr.sw <- data.frame(pop.code.st)
colnames(strata.hr.sw) <- "hr.sw"

## the actual data
sample.id <- read.gdsn(index.gdsn(wheat.data, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
# finding the informative subset of genotypes
snpset <- snpgdsSelectSNP(wheat.data, autosome.only = F, maf = 0.05, missing.rate = 0.1)
snpset.id <- unlist(snpset)
informative <- match(snpset.id, snp.id)
genotypes <- genotypes[informative,]
snpgdsClose(wheat.data)
# making the data palatable by genind
genotypes <- replace(genotypes, genotypes == 0, "A")
genotypes <- replace(genotypes, genotypes == 2, "B")
genotypes <- replace(genotypes, genotypes == 3, "N")

## making the genind objects
wheat.poppr.bp <- df2genind(t(data.frame(genotypes[,index.bp])), 
                         ind.names = as.character(sample.id)[index.bp], 
                         loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                         ncode = 1, strata = strata.bp)
wheat.poppr.or <- df2genind(t(data.frame(genotypes[,index.or])), 
                            ind.names = as.character(sample.id)[index.or], 
                            loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                            ncode = 1, strata = strata.or)
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
wheat.poppr.de <- df2genind(t(data.frame(genotypes[,index.de])), 
                            ind.names = as.character(sample.id)[index.de], 
                            loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                            ncode = 1, strata = strata.de)
wheat.poppr.mc <- df2genind(t(data.frame(genotypes[,index.mc])), 
                            ind.names = as.character(sample.id)[index.mc], 
                            loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                            ncode = 1, strata = strata.mc)
wheat.poppr.pk <- df2genind(t(data.frame(genotypes)), 
                            ind.names = as.character(sample.id), 
                            loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                            ncode = 1, strata = strata.pk)
wheat.poppr.hr.sw <- df2genind(t(data.frame(genotypes)), 
                            ind.names = as.character(sample.id), 
                            loc.names = snp.id[informative], NA.char = "N", ploidy = 1, type = "codom", 
                            ncode = 1, strata = strata.hr.sw)

# saving the genlight object
save(wheat.poppr.bp, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_bp.RData")
save(wheat.poppr.or, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_or.RData")
save(wheat.poppr.st, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_st.RData")
save(wheat.poppr.co, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_co.RData")
save(wheat.poppr.se, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_se.RData")
save(wheat.poppr.de, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_de.RData")
save(wheat.poppr.mc, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_mc.RData")
save(wheat.poppr.pk, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_pk.RData")
save(wheat.poppr.hr.sw, file = "C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genind_hr_sw.RData")
