library(tidyverse)
library(SNPRelate)
library(dbscan)
library(scrime)

## loading the gds of the data and pullling some attributes out
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds"
source("src\\R_functions\\data_loading.R")

genotypes <- genotypes %>%
    replace(. == 0, 1) %>%
    replace(. == 3, NA)

genotypes_imputed <- knncatimpute(t(genotypes))

full_hdbscan <- hdbscan(genotypes_imputed, minPts = 9)

write_rds(full_hdbscan,
    path = "Data\\Intermediate\\dbscan\\full_hdbscan.rds"
)

# examine group composition
# table(full_hdbscan$cluster)
# table(full_hdbscan$cluster[which(habit == "Winter")])
# table(desig[which(habit == "Winter")])
# table(full_hdbscan$cluster[which(desig == "HRS")])
# table(full_hdbscan$cluster[which(mc == "CWES")])
# table(desig[which(full_hdbscan$cluster == 0)])
# table(mc[which(full_hdbscan$cluster == 0)])
# table(bp[which(full_hdbscan$cluster == 0)])
# table(full_hdbscan$cluster[which(bp == "FOREIGN")])
# cbind(sample_id, as.character(desig), as.character(mc))[which(full_hdbscan$cluster == 0),]
# length(which(desig == "SWS"))