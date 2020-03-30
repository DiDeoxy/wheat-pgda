library(SNPRelate)

source("src/R_functions/funcs_gds_parse_create.R")

wheat_data <- parse_gds("phys_subset_sample")

wheat_gds <- snpgdsOpen("Data/Intermediate/GDS/phys_subset_sample.gds")

# ld pruned set of markes
set.seed(1000)
kept_id <- unlist(snpgdsLDpruning(wheat_gds,
  autosome.only = FALSE,
  maf = 0.05, slide.max.bp = 1e7, ld.threshold = 0.7
))
kept_index <- match(kept_id, wheat_data$snp$id)

snpgds_create_snp_subset(
  wheat_data, "phys_subset_sample_ld_pruned", kept_index
)