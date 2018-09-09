library(tidyverse)
library(SNPRelate)
library(dbscan)
library(scrime)

## loading the gds of the data and pullling some attributes out
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_both.gds"
source("src\\R_functions\\data_loading.R")

genotypes <- genotypes %>%
             replace(. == 0, 1) %>%
             replace(. == 3, NA)

genotypes_imputed <- knncatimpute(t(genotypes))

full_hdbscan <- hdbscan(genotypes_imputed, minPts = 9)

table(full_hdbscan$cluster)
# table(desig[which(full_hdbscan$cluster == 2)])
# table(desig)

write_rds(full_hdbscan,
          path = "Data\\Intermediate\\dbscan\\full_hdbscan.rds")