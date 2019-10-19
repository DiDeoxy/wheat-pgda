library(tidyverse)
library(SNPRelate)

source("src\\R_functions\\funcs_gds_parse_create.R")

cluster <- read_rds("Data\\Intermediate\\hdbscan\\wheat_hdbscan.rds")$cluster

subsets <- list(type = c("HRS", "HRW", "SWS"), cluster = c(5, 1, 2))

# creating gds subsets for major phenotypic/cluster types in the data set
# used in analyses of diversity within these groups
for (i in 1:length(subsets$type)) {
  wheat_data <- parse_gds("phys_subset_sample")
  clustered <- which(
    wheat_data$sample$annot$pheno == subsets$type[[i]] &
      cluster == subsets$cluster[[i]]
  )
  for (name in names(wheat_data$sample$annot)) {
    wheat_data$sample$annot[[name]] <- factor(
      wheat_data$sample$annot[[name]][clustered]
    )
  }
  snpgds_create_snp_subset(
    wheat_data, str_c("C", subsets$type[[i]], "_phys_subset_sample_pruned"),
    clustered
  )
}