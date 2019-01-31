library(tidyverse)
library(mmod)

locus_by_locus_comps <- c(
  "chrs_chrw_csws"
)[9]

m = 2/3
k = 3
km <- k * m
k_msqr <- k * m^2
km_i <- as.integer(km)
km_d <- (km - km_i)

max_D <- (2 * k * (km_i + (km_d)^2 - k_msqr)) / ((k - 1) * (k - 2 * km_d * (1 - km_d)))

for (comp in locus_by_locus_comps) {
  print(comp)
  comp_genind <- read_rds(
    "Data/Intermediate/Adegenet/chrs_chrw_csws_genind.rds"
  )
  write_rds(
    D_Jost(comp_genind)[[1]] / max_D,
    str_c("Data/Intermediate/mmod/", comp, "_Jost_D.rds")
  )
}