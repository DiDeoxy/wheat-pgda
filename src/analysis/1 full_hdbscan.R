library(SNPRelate)
library(plyr)
# install.packages("dbscan")
library(dbscan)
library(scrime)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the gds of the data and pullling some attributes out
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes <- replace(genotypes, genotypes == 3, NA)

genotypesImputed <- knncatimpute(t(genotypes))

wheatHdbscan <- hdbscan(genotypesImputed, minPts = 9)

save(wheatHdbscan, file = "Data\\Intermediate\\dbscan\\wheatHdbscan.RData")

# cluster <- wheatHdbscan$cluster
# table(cluster)
# table(mc[which(cluster == 3)])
# table(cluster[which(mc == "CPSR")])
