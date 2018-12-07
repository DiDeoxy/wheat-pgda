library(adegenet)

# find the underlying genetoypes of the extreme markers near genes

wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

################################################################################
# Lr10
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 1 &
      wheat_data$snp$pos_mb >= 8.296948 &
      wheat_data$snp$pos_mb <= 11.571978
  ),
] %>% nrow()

comp_genind <- read_rds(str_c("Data/Intermediate/Adegenet/Lr10_genind.rds"))

markers <- c(
  "RAC875_c57939_78"
)

for (marker in markers) {
  print(marker)
  tibble(
      pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

################################################################################
# Lr21
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 3 &
      wheat_data$snp$pos_mb >= 0.578085 &
      wheat_data$snp$pos_mb <= 2.497367
  ),
] %>% nrow()

comp_genind <- read_rds(str_c("Data/Intermediate/Adegenet/Lr21_genind.rds"))

markers <- c("D_contig09376_670")

for (marker in markers) {
  print(marker)
  tibble(
      pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

################################################################################
# lr22a
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 6 &
      wheat_data$snp$pos_mb >= 16.339182 &
      wheat_data$snp$pos_mb <= 16.356598
  ),
] %>% nrow()

comp_genind <- read_rds(str_c("Data/Intermediate/Adegenet/Lr22a_genind.rds"))

markers <- c(
  "wsnp_Ra_c25656_35227058"
)

for (marker in markers) {
  print(marker)
  tibble(
      pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

################################################################################
# lr1
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 15 &
      wheat_data$snp$pos_mb >= 561.705258 &
      wheat_data$snp$pos_mb <= 561.897587
  ),
] %>% nrow()

comp_genind <- read_rds(str_c("Data/Intermediate/Adegenet/Lr1_genind.rds"))

markers <- c(
  "wsnp_Ex_c11055_17927668"
)

for (marker in markers) {
  print(marker)
  tibble(
      pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}

################################################################################
# Lr34
wheat_data$snp[
  which(
    wheat_data$snp$chrom == 21 &
      wheat_data$snp$pos_mb >= 47.402696 &
      wheat_data$snp$pos_mb <= 58.644217
  ),
] %>% nrow()

comp_genind <- read_rds(str_c("Data/Intermediate/Adegenet/Lr34_genind.rds"))

markers <- c("Kukri_c32845_116")

for (marker in markers) {
  print(marker)
  tibble(
      pop = comp_genind$pop, allele = comp_genind$tab[, str_c(marker, ".A")]
  ) %>% table() %>% print()
}