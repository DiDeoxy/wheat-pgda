# import file paths and functions
source("wheat-pgda/src/R/file_paths.R")
import::from(dbscan, "hdbscan")
import::from(magrittr, "%>%")
import::from(pgda, "snpgds_parse")
import::from(readr, "write_csv", "write_rds")
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

################################################################################
# output some tables

write.csv(
  table(
    as.character(wheat_hdbscan$cluster),
    as.character(wheat_data$sample$annot$mtg)
  ) %>% as.data.frame.matrix(),
  "results/clusters_vs_mtg.csv"
)

write.csv(
  table(
    as.character(wheat_hdbscan$cluster),
    as.character(wheat_data$sample$annot$mc)
  ) %>% as.data.frame.matrix(),
  "results/clusters_vs_mc.csv"
)
