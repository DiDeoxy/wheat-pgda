library(tidyverse)
library(SNPRelate)
library(adegenet)
library(plyr)

source("src\\R_functions\\funcs_gds_parse_create.R")

# load the data
wheat_data <- parse_gds("phys_subset_sample_ld_pruned")
cluster <- read_rds("Data\\Intermediate\\hdbscan\\wheat_hdbscan.rds")$cluster

# making the genotype data palatable by genind
wheat_data$genotypes <- wheat_data$genotypes %>%
  replace(. == 0, "A") %>%
  replace(. == 2, "B") %>%
  replace(. == 3, "")

# make unknown and small groups uniformaly NA
wheat_data$sample$annot$bp <- revalue(
  wheat_data$sample$annot$bp,
  replace = (c("FOREIGN" = "N/A", "USA" = "N/A")
  )
)

# create vecotrs giving membership or not of varieties to each of the three
# largest phenotypic groups
major_phenos <- list()
for (pheno in c("HRS", "HRW", "SWS")) {
  major_phenos[[pheno]] <- as.character(wheat_data$sample$annot$pheno)
  subset_index <- which(wheat_data$sample$annot$pheno == pheno)
  major_phenos[[pheno]][-subset_index] <- "Other"
}

# create membership vectors of each cluster vs all other
cluster_isolated <- list()
for (cluster_num in c(1, 2, 3, 4, 5)) {
  cluster_isolated[[cluster_num]] <- cluster
  cluster_isolated[[cluster_num]][-which(cluster == cluster_num)] <- 0
}

# create a data frame of the different strata
strata <- tibble(
  era = wheat_data$sample$annot$era,
  bp = wheat_data$sample$annot$bp,
  texture = wheat_data$sample$annot$texture,
  colour = wheat_data$sample$annot$colour,
  habit = wheat_data$sample$annot$habit,
  pheno = wheat_data$sample$annot$pheno,
  hrs = as.factor(major_phenos$HRS),
  hrw = as.factor(major_phenos$HRW),
  sws = as.factor(major_phenos$SWS),
  clusters = as.factor(cluster),
  cluster1 = as.factor(cluster_isolated[[1]]),
  cluster2 = as.factor(cluster_isolated[[2]]),
  cluster3 = as.factor(cluster_isolated[[3]]),
  cluster4 = as.factor(cluster_isolated[[4]]),
  cluster5 = as.factor(cluster_isolated[[5]])
)

# turn the genotype data and strate into a genind
wheat_genind <- df2genind(
  t(wheat_data$genotypes),
  ind.names = wheat_data$sample$id,
  loc.names = wheat_data$snp$id, ploidy = 1, type = "codom", ncode = 1,
  strata = strata
)

write_rds(wheat_genind,
  path = paste0("Data\\Intermediate\\Adegenet\\wheat_genind.rds")
)