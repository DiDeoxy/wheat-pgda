library(tidyverse)

# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/funcs_locus_by_locus.R")

# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

# find the max position of any marker on each genome for xlims
max_genome_lengths <- calc_max_genome_lengths(wheat_data)

# names of groups to be plotted
groups <- c("chrs_csws", "chrs_chrw", "csws_chrw")

# find the jost'd values of each marker in each Gene and add to data set
wheat_data <- add_group_stat(wheat_data, groups)

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

write_csv(jost_d_stats[chroms_arranged, ], "Results/loci/D/D_stats.csv")