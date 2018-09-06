library(tidyverse)
library(SNPRelate)
library(dbscan)
library(scrime)

## loading the gds of the data and pullling some attributes out
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("src\\functions\\data_loading.R")

genotypes <- genotypes %>%
             replace(. == 0, 1) %>%
             genotypes(. == 3, NA)

genotypesImputed <- knncatimpute(t(genotypes))

wheatHdbscan <- hdbscan(genotypesImputed, minPts = 9)

save(wheatHdbscan, file = "Data\\Intermediate\\dbscan\\wheatHdbscan.RData")

# cluster <- wheatHdbscan$cluster
# table(cluster)
# table(mc[which(cluster == 3)])
# table(cluster[which(mc == "CPSR")])
