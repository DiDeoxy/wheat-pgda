## outputs data in a format readable by STRUCTURE
library(SNPRelate)
library(adegenet)
library(poppr)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))

# making the data palatable by genind
genotypes <- replace(genotypes, genotypes == 0, "A")
genotypes <- replace(genotypes, genotypes == 2, "B")
genotypes <- replace(genotypes, genotypes == 3, "N")

pop.code.bp <- read.gdsn(index.gdsn(wheat, "samp.annot/BP"))
snpgdsClose(wheat)

## for making subsets of the samples for the different categories excluding unknown samples
# index.bl <- which(pop.code.bp == "AAFC BEAVERLODGE RS")
index.cr <- which(pop.code.bp == "AAFC CEREAL RC")
# index.lc <- which(pop.code.bp == "AAFC LACOMBE RC")
index.lb <- which(pop.code.bp == "AAFC LETHBRIDGE RC")
# index.sc <- which(pop.code.bp == "AAFC SWIFT CURRENT")
# index.wi <- which(pop.code.bp == "AAFC WINNIPEG")
# index.au <- which(pop.code.bp == "AGRICORE UNITED/AGRIPRO")
# index.ap <- which(pop.code.bp == "AGRIPRO")
index.cdc <- which(pop.code.bp == "CROP DEVELOPMENT CENTRE")
# index.swp <- which(pop.code.bp == "SASKATCHEWAN WHEAT POOL")
index.spa <- which(pop.code.bp == "SPA RC")
index.sy <- which(pop.code.bp == "SYNGENTA")
# index.ugg <- which(pop.code.bp == "UNITED GRAIN GROWERS/AGRIPRO")
# index.ua <- which(pop.code.bp == "UNIVERSITY OF ALBERTA")
# index.um <- which(pop.code.bp == "UNIVERSITY OF MANITOBA")


## strata
strata.cr <- as.data.frame(rep(0, 357))
strata.cr[index.cr,] <- pop.code.bp[index.cr]
colnames(strata.cr) <- "cereal_none"

strata.lb <- as.data.frame(rep(0, 357))
strata.lb[index.lb,] <- pop.code.bp[index.lb]
colnames(strata.lb) <- "lethbridge_none"

strata.cdc <- as.data.frame(rep(0, 357))
strata.cdc[index.cdc,] <- pop.code.bp[index.cdc]
colnames(strata.cdc) <- "cdc_none"

strata.spa <- as.data.frame(rep(0, 357))
strata.spa[index.spa,] <- pop.code.bp[index.spa]
colnames(strata.spa) <- "spa_none"

strata.sy <- as.data.frame(rep(0, 357))
strata.sy[index.sy,] <- pop.code.bp[index.sy]
colnames(strata.sy) <- "syngenta_none"

## making the genind objects
wheat.poppr.cr <- df2genind(t(data.frame(genotypes)), 
                            ind.names = sample.id, loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.cr)
wheat.poppr.lb <- df2genind(t(data.frame(genotypes)), 
                            ind.names = sample.id, loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.lb)
wheat.poppr.cdc <- df2genind(t(data.frame(genotypes)), 
                            ind.names = sample.id, loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.cdc)
wheat.poppr.spa <- df2genind(t(data.frame(genotypes)), 
                            ind.names = sample.id, loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.spa)
wheat.poppr.sy <- df2genind(t(data.frame(genotypes)), 
                            ind.names = sample.id, loc.names = snp.id, NA.char = "N", 
                            ploidy = 1, type = "codom", ncode = 1, strata = strata.sy)

## amova
wheat.amova <- poppr.amova(wheat.poppr.cr, hier = ~cereal_none, missing = "genotype", within = F, clonecorrect = F)
wheat.amova

wheat.amova <- poppr.amova(wheat.poppr.lb, hier = ~lethbridge_none, missing = "genotype", within = F, clonecorrect = F)
wheat.amova

wheat.amova <- poppr.amova(wheat.poppr.cdc, hier = ~agripro_none, missing = "genotype", within = F, clonecorrect = F)
wheat.amova

wheat.amova <- poppr.amova(wheat.poppr.spa, hier = ~spa_none, missing = "genotype", within = F, clonecorrect = F)
wheat.amova

wheat.amova <- poppr.amova(wheat.poppr.sy, hier = ~syngenta_none, missing = "genotype", within = F, clonecorrect = F)
wheat.amova
