library(tidyverse)

# load custom functions
source("src\\R_functions\\funcs_gds_parse_create.R")

# load the data from the gds object
wheat_data <- parse_gds("phys_subset_sample")

# add columns of the D values of each comparison
for (group in c("chrs_csws", "chrs_chrw", "csws_chrw")) {
    comp_Jost_D <- read_rds(str_c(
        "Data\\Intermediate\\mmod\\", group,
        "_Jost_D.rds"
    ))[[1]]
    wheat_data$snp <- wheat_data$snp %>% add_column(!!group := comp_Jost_D)
}

chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()

jost_d_stats <- by(
  wheat_data$snp, wheat_data$snp$chrom, function(chrom) {
    chrom %>% summarise(
      chrs_csws_mean = mean(chrs_csws, na.rm = TRUE),
      chrs_csws_sd = sd(chrs_csws, na.rm = TRUE),
      chrs_chrw_mean = mean(chrs_chrw, na.rm = TRUE),
      chrs_chrw_sd = sd(chrs_chrw, na.rm = TRUE),
      csws_chrw_mean = mean(csws_chrw, na.rm = TRUE),
      csws_chrw_sd = sd(csws_chrw, na.rm = TRUE)
    )
  }
) %>%
  do.call(rbind, .) %>%
  add_row(
    chrs_csws_mean = mean(wheat_data$snp$chrs_csws),
    chrs_csws_sd = sd(wheat_data$snp$chrs_csws),
    chrs_chrw_mean = mean(wheat_data$snp$chrs_chrw),
    chrs_chrw_sd = sd(wheat_data$snp$chrs_chrw),
    csws_chrw_mean = mean(wheat_data$snp$csws_chrw),
    csws_chrw_sd = sd(wheat_data$snp$csws_chrw)
  ) %>%
  add_column(Chromosome = c(chroms, "All")) %>%
  select(Chromosome, everything())

chroms_arranged <- c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3), 22)

write_csv(jost_d_stats[chroms_arranged, ], "Results\\loci\\D\\D_stats.csv")