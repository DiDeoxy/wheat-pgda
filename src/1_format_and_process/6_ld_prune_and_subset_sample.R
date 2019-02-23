library(SNPRelate)
library(tidyverse)

source("src/R_functions/funcs_gds_parse_create.R")

# create ld pruned set of markes
wheat_gds <- snpgdsOpen(
  "Data/Intermediate/GDS/maf_and_mr_pruned_phys_sample_subset.gds"
)
set.seed(1000)
kept_id <- snpgdsLDpruning(
  wheat_gds, autosome.only = FALSE, slide.max.bp = 1e7, ld.threshold = 0.7
) %>% unlist()
snpgdsClose(wheat_gds)

wheat_data <- parse_gds("maf_and_mr_pruned_phys_sample_subset")
kept_index <- match(kept_id, wheat_data$snp$id)

snpgds_create_snp_subset(
  wheat_data, "ld_pruned_phys_sample_subset", kept_index
)