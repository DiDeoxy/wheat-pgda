library(tidyverse)
library(poppr)
# install.packages("mmod")
library(mmod)

comparisons <- c(
  "chrs_csws", "chrs_chrw", "csws_chrw", "Lr1", "Lr10", "Lr21", "Lr22a", "Lr34"
)

for (comparison in comparisons) {
  comp_genind <- read_rds(
    str_c("Data\\Intermediate\\Adegenet\\", comparison, "_genind.rds")
  )
  write_rds(
    diff_stats(comp_genind, phi_st = TRUE),
    str_c("Data\\Intermediate\\mmod\\", comparison, "_diff_stats.rds")
  )
}

comp_genind <- read_rds(
  str_c("Data\\Intermediate\\Adegenet\\Lr1_genind.rds")
)
comp_genind$strata
comp_genind$tab[1:10, 1:10]
head(read_rds("Data\\Intermediate\\mmod\\Lr1_diff_stats.rds")[[1]])