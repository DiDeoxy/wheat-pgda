# import file paths and functions
source(file.path("src", "R", "file_paths.R"))
import::from(dbscan, "hdbscan")
import::from(magrittr, "%>%")
import::from(pgda, "snpgds_parse")
import::from(readr, "write_rds")
import::from(scrime, "knncatimpute")

# load the ld pruned wheat data
wheat_data <- snpgds_parse(ld_phys_gds)

# impute the missing data
wheat_imputed <- wheat_data$genotypes %>%
    replace(. == 0, values = 1) %>%
    replace(. == 3, values = NA) %>%
    t() %>%
    knncatimpute()

# hdbscan the data
wheat_hdbscan <- wheat_imputed %>% hdbscan(minPts = 9)

# write the data out
write_rds(wheat_hdbscan, path = hdbscan)