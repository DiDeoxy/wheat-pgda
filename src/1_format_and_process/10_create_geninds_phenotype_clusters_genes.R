library(tidyverse)
library(plyr)
library(SNPRelate)
library(adegenet)
library(dbscan)

source("src\\R_functions\\funcs_gds_parse_create.R")

wheat_data <- parse_gds("phys_subset_sample")
cluster <- read_rds("Data\\Intermediate\\hdbscan\\wheat_hdbscan.rds")$cluster

# making the genotype data palatable by genind
wheat_data$genotypes <- wheat_data$genotypes %>%
  replace(. == 0, "A") %>%
  replace(. == 2, "B") %>%
  replace(. == 3, "")

# connect the phenotypes with their major clusters
subsets <- list(pheno = c("HRS", "HRW", "SWS"), cluster = c(5, 1, 2))

# making indices in order to create subsets containing just the varieties of
# the phenotype/cluster groups we are interested in
subset_indices <- list()
for (i in 1:length(subsets[[1]])) {
  subset_indices[[str_c("C", subsets$pheno[[i]])]] <- which(
    wheat_data$sample$annot$pheno == subsets$pheno[[i]] &
    cluster == subsets$cluster[[i]]
  )
}
index_chrs_csws <- c(subset_indices$CHRS, subset_indices$CSWS)
index_chrs_chrw <- c(subset_indices$CHRS, subset_indices$CHRW)
index_csws_chrw <- c(subset_indices$CSWS, subset_indices$CHRW)

grouping <- list(
  list(index_chrs_csws, index_chrs_chrw, index_csws_chrw),
  list("chrs_csws", "chrs_chrw", "csws_chrw")
)

for (i in 1:length(grouping[[1]])) {
  index <- grouping[[1]][[i]]
  name <- grouping[[2]][[i]]

  strata <- data.frame(as.character(wheat_data$sample$annot$pheno[index]))
  colnames(strata) <- name

  subset_genind <- df2genind(t(data.frame(wheat_data$genotypes[, index])),
    ind.names = wheat_data$sample$id[index], loc.names = wheat_data$snp$id,
    ploidy = 1, type = "codom", ncode = 1,
    strata = strata, pop = strata[[1]]
  )

  write_rds(subset_genind,
    path = str_c("Data\\Intermediate\\Adegenet\\", name, "_genind.rds")
  )
}

# making subsetted geninds of groups containing and not containing alleles of
# certain lR genes
gene_pres <- read_csv(
  "Data\\Intermediate\\Aligned_genes\\gene_presence_randhawa.csv",
  col_names = c("sample", "gene_A", "gene_B", "gene_C")
)
for (gene in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
  carriers <- c(which(gene_pres$gene_A == gene), 
    which(gene_pres$gene_B == gene), which(gene_pres$gene_C == gene)
  )
  index_carriers <- which(gene_pres$sample[carriers] %in% wheat_data$sample$id)
  index_not_carriers <- 
    which(gene_pres$sample[-carriers] %in% wheat_data$sample$id)
  num_carriers <- length(index_carriers)

  strata <- data.frame(
    c(
      rep(gene, length(index_carriers)),
      rep(str_c("not_", gene), length(index_not_carriers))
    )
  )
  colnames(strata) <- gene

  subset_genind <- df2genind(
    t(wheat_data$genotypes[,c(index_carriers, index_not_carriers)]),
    ind.names = 
      as.character(wheat_data$sample$id)[c(index_carriers,index_not_carriers)],
    loc.names = wheat_data$snp$id, ploidy = 1, type = "codom",
    ncode = 1, strata = strata, pop = strata[[1]]
  )

  write_rds(subset_genind,
    path = str_c("Data\\Intermediate\\Adegenet\\", gene, "_genind.rds"))
}
