library(tidyverse)
library(SNPRelate)

source("src\\R_functions\\funcs_gds_parse_create.R")

# load the data
wheat_data <- parse_gds("phys_subset_sample")

cluster <- read_rds("Data\\Intermediate\\dbscan\\full_hdbscan.rds")$cluster

# creating gds subsets for major phenotypic/cluster types in the data set
# used in analyses of diversity within these groups
hrs <- which(desig == "HRS")
chrs <- hrs[which(hrs %in% which(cluster == 5))]

chrs_sample_annot <- list()
for (name in names(wheat_data$sample$annot)) {
  chrs_sample_annot[[name]] <- factor(sample_annot[[name]][chrs])
}

snpgdsCreateGeno(
  "Data\\Intermediate\\GDS\\chrs_phys_subset_sample_pruned_floor.gds",
  genmat = genotypes[, chrs],
  sample.id = sample_id[chrs],
  snp.id = snp_id,
  snp.chromosome = snp_chrom,
  snp.position = snp_pos,
  other.vars = list(sample_annot = chrs_sample_annot),
  snpfirstdim = T
)

# all HRW varieites are in cluster 1
hrw <- which(desig == "HRW")
chrw <- hrs[which(hrw %in% which(cluster == 1))]
chrw_sample_annot <- list()
for (name in names(sample_annot)) {
  chrw_sample_annot[[name]] <- factor(sample_annot[[name]][chrw])
}
snpgdsCreateGeno(
  "Data\\Intermediate\\GDS\\chrw_phys_subset_sample_pruned_floor.gds",
  genmat = genotypes[, chrw],
  sample.id = sample_id[chrw],
  snp.id = snp_id,
  snp.chromosome = snp_chrom,
  snp.position = snp_pos,
  other.vars = list(sample_annot = chrw_sample_annot),
  snpfirstdim = T
)

# SWS spring varietes are spread over several clusters but the vast mjority
# are in cluster 2
sws <- which(desig == "SWS")
csws <- sws[which(sws %in% which(cluster == 3))]
csws_sample_annot <- list()
for (name in names(sample_annot)) {
  csws_sample_annot[[name]] <- factor(sample_annot[[name]][csws])
}
snpgdsCreateGeno(
  "Data\\Intermediate\\GDS\\csws_phys_subset_sample_pruned_floor.gds",
  genmat = genotypes[, csws],
  sample.id = sample_id[csws],
  snp.id = snp_id,
  snp.chromosome = snp_chrom,
  snp.position = snp_pos,
  other.vars = list(sample_annot = csws_sample_annot),
  snpfirstdim = T
)