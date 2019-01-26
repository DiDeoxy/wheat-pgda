library(tidyverse)
library(SNPRelate)

source("src\\R_functions\\funcs_gds_parse_create.R")

wheat_data <- parse_gds("phys_subset_sample")

wheat_data$genotypes <- wheat_data$genotypes %>%
  replace(. == 0, 1) %>%
  replace(. == 3, 0)

gene_pres <- read_csv(
  "Data\\Intermediate\\Aligned_genes\\gene_presence_randhawa.csv",
  col_names = c("sample", "gene_A", "gene_B", "gene_C")
)

make_table <- function (id, geno, indiv, num_carriers) {
  if (index <- num_carriers) {
    c(individual = str_c("P1.", indiv), population = "P1",
      fragment.length = geno, locus = id)
  } else {
    c(individual = str_c("P2.", indiv - num_carriers),
      population = "P2", fragment.length = snp_geno[indiv], locus = id)
  }
}

make_table_caller <- function (snp, num_carriers) {
  snp_id <- snp[1]
  snp_geno <- snp[2:length(snp)]
  snp_table <- tibble(individual = character(), population = character(),
    fragment.length = integer(), locus = character())
  mapply(make_table, rep(snp_id, length(snp_geno)), snp_geno,
    1:length(snp_geno), rep(num_carriers, length(snp_geno)), USE.NAMES = FALSE
  )
}

for (gene in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
  carriers <- c(which(gene_pres$gene_A == gene), 
    which(gene_pres$gene_B == gene), which(gene_pres$gene_C == gene)
  )
  index_carriers <- which(gene_pres$sample[carriers] %in% wheat_data$sample$id)
  index_not_carriers <- 
    which(gene_pres$sample[-carriers] %in% wheat_data$sample$id)
  num_carriers <- length(index_carriers)

  comp_table <- 
    apply(
      cbind(wheat_data$snp$id, 
        wheat_data$genotypes[, c(index_carriers, index_not_carriers)]),
      1, num_carriers = num_carriers, make_table_caller
    ) %>% matrix(
      nrow = length(c(index_carriers, index_not_carriers)) * 
        nrow(wheat_data$genotypes),
      byrow = TRUE
    ) %>% as.tibble()
  colnames(comp_table) <- 
    c("individual", "population", "fragment.length", "locus")
  comp_table <- comp_table[-which(comp_table$fragment.length == "0"), ]
  
  write_tsv(comp_table, 
    str_c("Data\\Intermediate\\DEMEtics\\", gene, "_genind.tsv"),
  )
}
