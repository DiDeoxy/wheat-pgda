# import file paths and functions
source(file.path("src", "R", "file_paths.R"))
library(tidyverse)
library(pgda)

# mean_dist <- 14547261565 / 14003


coverage_by_chrom <- function (snp_data) {
  by(snp_data, snp_data$chrom, function (chrom_data) {
    intervals <- tibble(start = double(), end = double())
    cur_interval <- list(start = double(), end = double())
    for (i in 1:nrow(chrom_data)) {
      if (i == 1) {
        cur_interval$start <- max(0, chrom_data[i, ]$pos - half_mean_dist)
        cur_interval$end <- chrom_data[i, ]$pos + half_mean_dist
      } else if (i > 1 && i < nrow(chrom_data)) {
        if ((chrom_data[i, ]$pos - half_mean_dist) <= cur_interval$end) {
          cur_interval$end <- chrom_data[i, ]$pos + half_mean_dist
        } else {
          intervals <- intervals %>%
            add_row(start = cur_interval$start, end = cur_interval$end)
          cur_interval$start <- chrom_data[i, ]$pos - half_mean_dist
          cur_interval$end <- chrom_data[i, ]$pos + half_mean_dist
        }
      } else {
        if ((chrom_data[i, ]$pos - half_mean_dist) <= cur_interval$end) {
          cur_interval$end <- chrom_data[i, ]$pos + half_mean_dist
          intervals <- intervals %>%
            add_row(start = cur_interval$start, end = cur_interval$end)
        } else {
          intervals <- intervals %>%
            add_row(start = cur_interval$start, end = cur_interval$end)
          cur_interval$start <- chrom_data[i, ]$pos - half_mean_dist
          cur_interval$end <- chrom_data[i, ]$pos + half_mean_dist
          intervals <- intervals %>%
            add_row(start = cur_interval$start, end = cur_interval$end)
        }
      }
    }
    (intervals[, 2] - intervals[, 1]) %>% sum()
  }) %>% as.list() %>% unlist()
}

wheat_data <- snpgds_parse(file.path(gds, "maf_mr_filtered_phys.gds"))
half_mean_dist <- 0.5 * (wheat_data$chrom_lengths %>% sum() / 14003)
cov_by_chrom <- coverage_by_chrom(wheat_data$snp)
cov_by_chrom / wheat_data$chrom_lengths
cov_by_chrom %>% sum() / wheat_data$chrom_lengths %>% sum()

(wheat_data$chrom[seq(3, 21, 3)] %>% sum() / 1e6) /
((wheat_data$chrom[seq(1, 19, 3)] %>% sum() / 1e6) +
(wheat_data$chrom[seq(2, 20, 3)] %>% sum() / 1e6) +
(wheat_data$chrom[seq(3, 21, 3)] %>% sum() / 1e6))

(wheat_data$chrom[seq(1, 19, 3)] %>% sum() / 1e6) / 5567
(wheat_data$chrom[seq(2, 20, 3)] %>% sum() / 1e6) / 7266
(wheat_data$chrom[seq(3, 21, 3)] %>% sum() / 1e6) / 1170

wheat_data <- snpgds_parse(file.path(gds, "maf_mr_filtered_gen.gds"))
half_mean_dist <- 0.5 * (by(wheat_data$snp$pos, wheat_data$snp$chrom, max) %>% as.list() %>% unlist() / 14003)
cov_by_chrom <- coverage_by_chrom(wheat_data$snp)
cov_by_chrom / (by(wheat_data$snp$pos, wheat_data$snp$chrom, max) %>% as.list() %>% unlist())
cov_by_chrom %>% sum() / (by(wheat_data$snp$pos, wheat_data$snp$chrom, max) %>% as.list() %>% unlist()) %>% sum()


median(order_diffs)
getmode(order_diffs)
mean(order_diffs)
sd(order_diffs)
hist(order_diffs)
summary(order_diffs)

test1 <- c("a", "b", "c", "d", "e")
test2 <- c("c", "e", "b", "a", "d")

match(test1, test2)
match(rev(test1), rev(test2))
test2[match(test1, test2)]
rev(test2)[match(rev(test1), rev(test2))] %>% rev()