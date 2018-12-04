library(tidyverse)
library(SNPRelate)
library(dbscan)
library(scrime)

source("src/R_functions/funcs_gds_parse_create.R")

wheat_data <- parse_gds("ld_pruned_phys_sample_subset")

wheat_imputed <- wheat_data$genotypes %>%
    replace(. == 0, 1) %>%
    replace(. == 3, NA) %>%
    t() %>%
    knncatimpute()

wheat_hdbscan <- wheat_imputed %>% hdbscan(minPts = 9)

write_rds(wheat_hdbscan,
    path = "Data/Intermediate/hdbscan/wheat_hdbscan.rds"
)

wheat_hdbscan <- read_rds("Data/Intermediate/hdbscan/wheat_hdbscan.rds")

clusters <- wheat_hdbscan$cluster %>%
  replace(. == 0, "Noise") 

table(data.frame(wheat_data$sample$annot$pheno, clusters)) %>%
  as.data.frame.matrix() %>%
  write.csv("Results/pheno_cluster.csv", quote = FALSE)
table(data.frame(wheat_data$sample$annot$mc, clusters)) %>%
  as.data.frame.matrix() %>%
  write.csv("Results/mc_cluster.csv", quote = FALSE)