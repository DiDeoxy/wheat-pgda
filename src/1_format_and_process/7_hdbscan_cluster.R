library(tidyverse)
library(SNPRelate)
library(dbscan)
library(scrime)

source("src/R_functions/funcs_gds_parse_create.R")

wheat_data <- parse_gds("phys_subset_sample")

wheat_imputed <- wheat_data$genotypes %>%
    replace(. == 0, 1) %>%
    replace(. == 3, NA) %>%
    t() %>%
    knncatimpute()

wheat_hdbscan <- wheat_imputed %>% hdbscan(minPts = 8)

write_rds(wheat_hdbscan,
    path = "Data/Intermediate/hdbscan/wheat_hdbscan.rds"
)

# examine group composition
table(wheat_hdbscan$cluster)

# cluster 1
table(wheat_data$sample$annot$pheno[wheat_hdbscan$cluster == 1])
table(wheat_hdbscan$cluster[wheat_data$sample$annot$pheno == "HRW"])
table(wheat_hdbscan$cluster[wheat_data$sample$annot$pheno == "HWW"])
table(wheat_hdbscan$cluster[wheat_data$sample$annot$pheno == "SRW"])
table(wheat_hdbscan$cluster[wheat_data$sample$annot$pheno == "SWW"])
table(wheat_hdbscan$cluster[wheat_data$sample$annot$habit == "Winter"])

# cluster 2
table(wheat_data$sample$annot$pheno[wheat_hdbscan$cluster == 2])
table(wheat_hdbscan$cluster[wheat_data$sample$annot$pheno == "SWS"])

# cluster 3
table(wheat_data$sample$annot$pheno[wheat_hdbscan$cluster == 3])
table(wheat_data$sample$annot$mc[wheat_hdbscan$cluster == 3])
table(wheat_data$sample$annot$mc[wheat_data$sample$annot$pheno == "HWS" &
  wheat_hdbscan$cluster == 3])

# cluster 4
table(wheat_data$sample$annot$pheno[wheat_hdbscan$cluster == 4])
table(wheat_data$sample$annot$mc[wheat_hdbscan$cluster == 4])
table(wheat_data$sample$annot$pheno[wheat_data$sample$annot$mc == "CWGP" &
  wheat_hdbscan$cluster == 4])
table(wheat_data$sample$annot$pheno[wheat_data$sample$annot$mc == "N/A" &
  wheat_hdbscan$cluster == 4])

# cluster 5
table(wheat_data$sample$annot$pheno[wheat_hdbscan$cluster == 5])

# Noise
table(wheat_data$sample$annot$pheno[wheat_hdbscan$cluster == 0])
table(wheat_data$sample$annot$mc[wheat_hdbscan$cluster == 0])
table(wheat_hdbscan$cluster[wheat_data$sample$annot$pheno == "HRS"])