library(adegenet)
library(tidyverse)

comp_genind <- read_rds(str_c(
  "Data\\Intermediate\\Adegenet\\Lr34_genind.rds"
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


comp_genind <- read_rds(str_c(
  "Data\\Intermediate\\Adegenet\\csws_chrw_genind.rds"
))

markers <- c("TA001775-0903")

for (marker in markers) {
  print(marker)
  cbind(
    allele = comp_genind$tab[, str_c(marker, ".A")], comp_genind$strata
  ) %>% table() %>% print()
}