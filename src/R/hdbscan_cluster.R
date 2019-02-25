# load the ld pruned wheat data
wheat_data <- parse_gds(ld_gds)

# impute the missing data
wheat_imputed <- wheat_data$genotypes %>%
    replace(. == 0, 1) %>%
    replace(. == 3, NA) %>%
    t() %>%
    knncatimpute()

# hdbscan the data
wheat_hdbscan <- wheat_imputed %>% hdbscan(minPts = 9)

# write the data out
write_rds(wheat_hdbscan, path = hdbscan)