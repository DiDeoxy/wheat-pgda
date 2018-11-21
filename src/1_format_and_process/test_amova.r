library(adegenet)
library(tidyverse)

source("src/R_functions/funcs_gds_parse_create.R")

comp_genind <- read_rds(str_c(
  "Data/Intermediate/Adegenet/Lr34_genind.rds"
))

markers <- c(
  "Kukri_c92151_216", "Kukri_c32845_116",
  "TA002473-0717"
)

for (marker in markers) {
  print(marker)
  cbind(
    allele = comp_genind$tab[, str_c(marker, ".A")], comp_genind$strata
  ) %>% table() %>% print()
}




wheat_data <- parse_gds("phys_subset_sample")

wheat_data$snp[
  which(
    wheat_data$snp$chrom == 9 &
      wheat_data$snp$pos_mb >= 565.80153 &
      wheat_data$snp$pos_mb <= 576.550497
  ),
] %>% nrow()

which(wheat_data$snp$id == "wsnp_CAP7_c44_26549")
wheat_data$snp$pos[14260:14272]

comp_genind <- read_rds(str_c(
  "Data/Intermediate/Adegenet/csws_chrw_genind.rds"
))

markers <- c("BS00022055_51", "RAC875_c34939_963")

for (marker in markers) {
  print(marker)
  cbind(
    allele = comp_genind$tab[, str_c(marker, ".A")], comp_genind$strata
  ) %>% table() %>% print()
}

# CHRS
59 + 12 + 11 + 20 + 14 + 25 + 5 + 11 + 2 + 28 + 19 + 9
# CHRW
65 + 12 + 20 + 14 + 5 + 22 + 54 + 19
# CSWS
59 + 65 + 11 + 25 + 5 + 5 + 22 + 11 + 2 + 54 + 28 + 9
