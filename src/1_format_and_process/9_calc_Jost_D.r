library(tidyverse)
library(mmod)

m = 2/3
k = 3
km <- k * m
k_msqr <- k * m^2
km_i <- as.integer(km)
km_d <- (km - km_i)

max_D <- (2 * k * (km_i + (km_d)^2 - k_msqr)) / ((k - 1) * (k - 2 * km_d * (1 - km_d)))

comp <- "chrs_chrw_csws"

comp_genind <- read_rds(
  str_c("Data/Intermediate/Adegenet/", comp, "_genind.rds")
)
write_rds(
  D_Jost(comp_genind)[[1]] / max_D,
  str_c("Data/Intermediate/mmod/", comp, "_Jost_D.rds")
)