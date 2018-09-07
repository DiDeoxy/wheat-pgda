library(tidyverse)
library(SNPRelate)

# load the data
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\functions\\data_loading.R")

cluster <- read_rds("Data\\Intermediate\\dbscan\\full_hdbscan.rds")$cluster

# creating gds subsets for major phenotypic/cluster types in the data set
# used in analyses of diversity within these groups
hrs <- which(desig == "HRS")
chrs <- hrs[which(hrs %in% which(cluster == 5))]

chrs_samp_annot <- list()
for(name in names(samp_annot)) {
  chrs_samp_annot[[name]] <- factor(samp_annot[[name]][chrs])
}

snpgdsCreateGeno("Data\\Intermediate\\GDS\\chrs_phys_subset_sample.gds",
                 genmat = genotypes[,chrs],
                 sample.id = sample_id[chrs],
                 snp.id = snp_id,
                 snp.chromosome = snp_chrom,
                 snp.position = snp_pos,
                 other.vars = list(samp_annot = chrs_samp_annot),
                 snpfirstdim = T)

# all HRW varieites are in cluster 1
hrw <- which(desig == "HRW")
hrw_samp_annot <- list()
for(name in names(samp_annot)) {
  hrw_samp_annot[[name]] <- factor(samp_annot[[name]][hrw])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\hrw_phys_subset_sample.gds",
                 genmat = genotypes[,hrw],
                 sample.id = sample_id[hrw],
                 snp.id = snp_id,
                 snp.chromosome = snp_chrom,
                 snp.position = snp_pos,
                 other.vars = list(samp_annot = hrw_samp_annot),
                 snpfirstdim = T)

# SWS spring varietes are spread over several clusters but the vast mjority
# are in cluster 2
sws <- which(desig == "SWS")
csws <- sws[which(sws %in% which(cluster == 2))]
csws_samp_annot <- list()
for(name in names(samp_annot)) {
  csws_samp_annot[[name]] <- factor(samp_annot[[name]][csws])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\csws_phys_subset_sample.gds",
                 genmat = genotypes[,csws],
                 sample.id = sample_id[csws],
                 snp.id = snp_id,
                 snp.chromosome = snp_chrom,
                 snp.position = snp_pos,
                 other.vars = list(samp_annot = csws_samp_annot),
                 snpfirstdim = T)