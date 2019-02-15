library(SNPRelate)
library(tidyverse)

source("src/R_functions/funcs_gds_parse_create.R")

# find markers with a a missing rate below 0.1
wheat_gds <- snpgdsOpen("Data/Intermediate/GDS/full_phys.gds")
kept_id <- unlist(snpgdsSelectSNP(
  wheat_gds, autosome.only = F, missing.rate = 0.1
))
snpgdsClose(wheat_gds)

# reomve these markers from the phys map
wheat_data <- parse_gds("full_phys")
kept_index <- match(kept_id, wheat_data$snp$id)
snpgds_create_snp_subset(wheat_data, "mr_pruned_phys", kept_index)

# remove these markers from the gen map
wheat_data <- parse_gds("full_gen")
kept_index <- match(kept_id, wheat_data$snp$id) %>% sort()
snpgds_create_snp_subset(wheat_data, "mr_pruned_gen", kept_index)

# eliminate those individuals that show identity by state
# (IBS, fractional identity) greater than 0.99
wheat_gds <- snpgdsOpen(
  "Data/Intermediate/GDS/mr_pruned_phys.gds"
)
IBS <- snpgdsIBS(wheat_gds, autosome.only = F)
snpgdsClose(wheat_gds)

pairs <- which(IBS$ibs >= 0.99, arr.ind = T)

indices <- vector()
for (i in 1:dim(pairs)[1]) {
  if (pairs[i, 1] == pairs[i, 2]) {
    indices <- c(indices, i)
  }
}
pairs <- pairs[-indices, ]
pairs

# i take the pairs and use a perl function to find the connected graphs
# i.e. those sets of pairs that form a connected triangle, square or more

# eliminate individuals so that only one from each pair, triangle, or square
# remains
NILs <- c(
  "PT434", "BW811", "AC Minto 1", "Avocet 1", "BW275 1", "BW395",
  "PT616", "BW427 1", "BW492", "BW948", "Carberry 1", "CDC Stanley 1",
  "PT754", "SWS349", "Somerset 1", "Stettler 1", "AC Reed 1", "SWS87",
  "SWS241", "SWS345", "SWS363", "SWS390", "SWS408", "SWS410"
)

# find the indices of the NILs
sample_index <- match(NILs, wheat_data$sample$id)

# mr pruned phys map
wheat_data <- parse_gds("mr_pruned_phys")
# create gds object without the NILs
snpgds_create_sample_subset(
  wheat_data, "mr_pruned_phys_sample_subset", sample_index
)

# mr pruned gen map
wheat_data <- parse_gds("mr_pruned_gen")
snpgds_create_sample_subset(
  wheat_data, "mr_pruned_gen_sample_subset", sample_index
)

# identify snps with a maf below 0.05
wheat_gds <- snpgdsOpen(
  "Data/Intermediate/GDS/mr_pruned_phys_sample_subset.gds"
)
kept_id <- snpgdsSelectSNP(wheat_gds, maf = 0.05, autosome.only = F)
snpgdsClose(wheat_gds)

# reomve these markers from the phys map
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")
kept_index <- match(kept_id, wheat_data$snp$id)
snpgds_create_snp_subset(
  wheat_data, "maf_and_mr_pruned_phys_sample_subset", kept_index
)

# remove these markers from the gen map
wheat_data <- parse_gds("mr_pruned_gen_sample_subset")
kept_index <- match(kept_id, wheat_data$snp$id) %>% sort()
snpgds_create_snp_subset(
  wheat_data, "maf_and_mr_pruned_gen_sample_subset", kept_index
)