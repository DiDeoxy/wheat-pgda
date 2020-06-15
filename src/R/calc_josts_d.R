source(file.path("/workspace/repos/wheat-pgda/src/R/file_paths.R"))
import::from(mmod, "D_Jost")
import::from(readr, "read_rds", "write_rds")

# calc maximal josts D from Alcala & Rosenberg, 2018 so we can normalize
k = 3 # num groups
m = (k - 1) / k
km <- k * m
k_msqr <- k * m^2
km_i <- as.integer(km)
km_d <- (km - km_i)
max_D <- (2 * k * (km_i + (km_d)^2 - k_msqr)) /
    ((k - 1) * (k - 2 * km_d * (1 - km_d)))

# load the genind object
comp_genind <- read_rds(file.path(geninds, "chrs_chrw_csws.rds"))

norm_josts_d <- D_Jost(comp_genind)[[1]] / max_D
norm_josts_d[which(norm_josts_d < 0)] <- 0

# calculate normalized josts d values and wrtie the out
write_rds(norm_josts_d, josts_d)