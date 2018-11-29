library(tidyverse)
library(SNPRelate)
library(adegenet)
library(plyr)

source("src/R_functions/funcs_gds_parse_create.R")

################################################################################
# Make geninds for different groupings of the sample
wheat_data <- parse_gds("phys_subset_sample_ld_pruned")

cluster <- read_rds("Data/Intermediate/hdbscan/wheat_hdbscan.rds")$cluster

# making the genotype data palatable by genind
wheat_data$genotypes <- wheat_data$genotypes %>%
  replace(. == 0, "A") %>%
  replace(. == 2, "B") %>%
  replace(. == 3, "")

# make unknown and small groups uniformaly NA
wheat_data$sample$annot$bp <- revalue(
  wheat_data$sample$annot$bp,
  replace = c("FOREIGN" = "N/A", "USA" = "N/A")
)

# create vectors giving membership or not of varieties to each of the three
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

df2genind(
  t(wheat_data$genotypes),
  ind.names = wheat_data$sample$id,
  loc.names = wheat_data$snp$id, ploidy = 1, type = "codom", ncode = 1,
  strata = strata
) %>%
write_rds(path = str_c("Data/Intermediate/Adegenet/strata_genind.rds"))

# ################################################################################
# wheat_data <- parse_gds("phys_subset_sample")

# # making subsetted geninds of groups containing and not containing alleles of
# # certain lR genes
# gene_pres <- read_csv(
#   "Data/Intermediate/Aligned_genes/gene_presence_randhawa.csv",
#   col_names = c("sample", "gene_A", "gene_B", "gene_C")
# )

# for (gene in c("Lr10", "Lr21", "Lr22a", "Lr1", "Lr34")[2]) {
#   carriers <- c(which(gene_pres$gene_A == gene), 
#     which(gene_pres$gene_B == gene), which(gene_pres$gene_C == gene)
#   )
#   index_carriers <- which(gene_pres$sample[carriers] %in% wheat_data$sample$id)
  
#   index_not_carriers <-
#     which(gene_pres$sample[-carriers] %in% wheat_data$sample$id)
#   num_carriers <- length(index_carriers)

#   # print out how many of each kind there are for each gene
#   print(gene)
#   print(num_carriers)
#   print(length(index_not_carriers))

#   pop <- c(
#     rep(gene, length(index_carriers)),
#     rep(str_c("not_", gene), length(index_not_carriers))
#   )

#   subset_genind <- df2genind(
#     t(wheat_data$genotypes[, c(index_carriers, index_not_carriers)]),
#     ind.names =
#       as.character(wheat_data$sample$id)[
#         c(index_carriers, index_not_carriers)
#       ],
#     loc.names = wheat_data$snp$id, ploidy = 1, type = "codom",
#     ncode = 1, pop = pop
#   )

#   write_rds(subset_genind,
#     path = str_c("Data/Intermediate/Adegenet/", gene, "_genind.rds"))
# }

# ################################################################################
# # making subset geninds of each comparison of the clustered varieties of the
# # three largest phenotypic groups
# wheat_data <- parse_gds("phys_subset_sample")

# # making indices in order to create subsets containing just the varieties of
# # the phenotype/cluster groups we are interested in
# pheno_indices <- list()
# pheno_indices[["CHRW"]] <- which(
#   wheat_data$sample$annot$pheno == "HRW" & cluster == 1
# )
# pheno_indices[["CSWS"]] <- which(
#   wheat_data$sample$annot$pheno == "SWS" & cluster == 2
# )
# pheno_indices[["CHRS"]] <- which(
#   wheat_data$sample$annot$pheno == "HRS" & cluster == 5
# )
# index_chrs_csws <- c(pheno_indices$CHRS, pheno_indices$CSWS)
# index_chrs_chrw <- c(pheno_indices$CHRS, pheno_indices$CHRW)
# index_csws_chrw <- c(pheno_indices$CSWS, pheno_indices$CHRW)

# grouping <- list(
#   list(index_chrs_csws, index_chrs_chrw, index_csws_chrw),
#   list("chrs_csws", "chrs_chrw", "csws_chrw")
# )

# for (i in 1:length(grouping[[1]])) {
#   index <- grouping[[1]][[i]]
#   name <- grouping[[2]][[i]]

#   pop <- as.character(wheat_data$sample$annot$pheno[index])

#   subset_genind <- df2genind(t(data.frame(wheat_data$genotypes[, index])),
#     ind.names = wheat_data$sample$id[index], loc.names = wheat_data$snp$id,
#     ploidy = 1, type = "codom", ncode = 1,
#     pop = pop
#   )

#   write_rds(subset_genind,
#     path = str_c("Data/Intermediate/Adegenet/", name, "_genind.rds")
#   )
# }