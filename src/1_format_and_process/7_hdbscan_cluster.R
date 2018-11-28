library(tidyverse)
library(SNPRelate)
library(dbscan)
library(scrime)

source("src/R_functions/funcs_gds_parse_create.R")

wheat_data <- parse_gds("phys_subset_sample_ld_pruned")
# length(wheat_data$snp$id)

wheat_imputed <- wheat_data$genotypes %>%
    replace(. == 0, 1) %>%
    replace(. == 3, NA) %>%
    t() %>%
    knncatimpute()

wheat_hdbscan <- wheat_imputed %>% hdbscan(minPts = 8)

write_rds(wheat_hdbscan,
    path = "Data/Intermediate/hdbscan/wheat_hdbscan.rds"
)

wheat_hdbscan <- read_rds("Data/Intermediate/hdbscan/wheat_hdbscan.rds")

clusters <- wheat_hdbscan$cluster %>%
  replace(. == 0, "Noise") 
  # %>%
  # replace(. == 1, "Cluster 1") %>%
  # replace(. == 2, "Cluster 2") %>%
  # replace(. == 3, "Cluster 3") %>%
  # replace(. == 4, "Cluster 4") %>%
  # replace(. == 5, "Cluster 5")

table(data.frame(wheat_data$sample$annot$pheno, clusters)) %>%
  as.data.frame.matrix() %>%
  write.csv("Results/pheno_cluster.csv", quote = FALSE)
table(data.frame(wheat_data$sample$annot$mc, clusters)) %>%
  as.data.frame.matrix() %>%
  write.csv("Results/mc_cluster.csv", quote = FALSE)