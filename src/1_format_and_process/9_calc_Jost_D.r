library(tidyverse)
library(mmod)

locus_by_locus_comps <- c(
  "Lr10", "Lr21", "Lr22a", "Lr1", "Lr34", "chrs_csws", "chrs_chrw", "csws_chrw",
  "chrs_chrw_csws"
)[9]

for (comp in locus_by_locus_comps) {
  print(comp)
  comp_genind <- read_rds(
    str_c("Data/Intermediate/Adegenet/", comp, "_genind.rds")
  )
  write_rds(
    D_Jost(comp_genind),
    str_c("Data/Intermediate/mmod/", comp, "_Jost_D.rds")
  )
}
pairwise_D(read_rds("Data/Intermediate/Adegenet/strata_genind.rds"))