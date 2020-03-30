# ################################################################################
wheat_data <- parse_gds(ld_gds)

# making the genotype data palatable by genind
wheat_data$genotypes <- wheat_data$genotypes %>%
  replace(. == 0, "A") %>%
  replace(. == 2, "B") %>%
  replace(. == 3, "")

# making subsetted geninds of groups containing and not containing alleles of
# certain lR genes
gene_pres <- read_csv(
  file.path(genes, "gene_presence_randhawa.csv"),
  col_names = c("sample", "gene_A", "gene_B", "gene_C")
)

for (gene in c("Lr10", "Lr21", "Lr22a", "Lr1", "Lr34")) {
  carriers <- c(which(gene_pres$gene_A == gene), 
    which(gene_pres$gene_B == gene), which(gene_pres$gene_C == gene)
  )
  index_carriers <- which(wheat_data$sample$id %in% gene_pres$sample[carriers])
  
  index_not_carriers <-
    which(wheat_data$sample$id %in% gene_pres$sample[-carriers])
  
  num_carriers <- length(index_carriers)

  # print out how many of each kind there are for each gene
  print(gene)
  print(num_carriers)
  print(length(index_not_carriers))

  pop <- c(
    rep(gene, length(index_carriers)),
    rep(str_c("not_", gene), length(index_not_carriers))
  )

  subset_genind <- df2genind(
    t(wheat_data$genotypes[, c(index_carriers, index_not_carriers)]),
    ind.names =
      as.character(wheat_data$sample$id)[
        c(index_carriers, index_not_carriers)
      ],
    loc.names = wheat_data$snp$id, ploidy = 1, type = "codom",
    ncode = 1, pop = pop
  )

  write_rds(subset_genind, path = file.path(geninds, str_c(gene, ".rds")))
}