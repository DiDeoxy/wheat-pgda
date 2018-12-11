library(tidyverse)

calc_eh <- function (genotypes) {
  apply(genotypes, 1, function (snp) {
    2 * ((sum(snp == 0) / sum(snp == 0 | 2)) *
      (sum(snp == 2) / sum(snp == 0 | 2)))
  })
}

calc_max_genome_lengths <- function(wheat_data) {
  by(wheat_data$snp$pos_mb, wheat_data$snp$chrom, max) %>% 
    (function (max_chrom_lengths) {
      c(
        A = max(max_chrom_lengths[seq(1, 19, 3)]),
        B = max(max_chrom_lengths[seq(2, 20, 3)]),
        D = max(max_chrom_lengths[seq(3, 21, 3)])
      )
    })
}

add_group_stat <- function(wheat_data, groups) {
  for (group in groups) {
    wheat_data$snp <- wheat_data$snp %>% 
      add_column(
        !!group := read_rds(
          str_c("Data/Intermediate/mmod/", group,"_Jost_D.rds")
        )[[1]]
      )
  }
  wheat_data
}
group <- "Lr34"

calc_extremes <- function(wheat_data, groups, prob = 0.975) {
  temp <- vector()
  for (group in groups) {
    temp[group] <- quantile(wheat_data$snp[[group]], prob = prob, na.rm = TRUE)
  }
  temp
}

calc_group_extreme_freqs <- function(wheat_data, extremes, prune = FALSE) {
  by(wheat_data$snp, wheat_data$snp$chrom,
    function(snp_data) {
      # initialize an empty tibble with columns for the important data
      ret <- tibble(
        chrom = integer(), group = character(), pos_mb = double(),
        freq = double(), extreme_D = double(), id = character()
      )
      # 
      for (group in names(extremes)) {
        for (i in 1:length(snp_data$id)) {
        num_nearby <- which(
          snp_data$pos_mb >= snp_data$pos_mb[i] - 5 &
          snp_data$pos_mb <= snp_data$pos_mb[i] + 5
        ) %>% length()
        freq <- integer()
        if (num_nearby >= 5 && i >= 3 && i <= (length(snp_data$id) - 2)) {
          freq <- sum(
            snp_data[[group]][(i - 2):(i + 2)] > extremes[[group]]
          ) / 5
        } else if (
          (i < 3 || i > (length(snp_data$id) - 2)) &&
          snp_data[[group]][i] > extremes[[group]]
        ) {
          freq <- 0.1
        } else {
          freq <- NA
        }
        if (! is.na(freq) && freq == 0) {
          freq <- NA
        }
          ret <- ret %>%
            add_row(
              chrom = snp_data$chrom[i], group = group,
              pos_mb = snp_data$pos_mb[i], freq = freq,
              extreme_D = ifelse(
                snp_data[[group]][i] > extremes[[group]], snp_data[[group]][i],
                NA
              ),
              id = snp_data$id[i]
            )
        }
      }
      # set all markers that don't have the highest frequency of extreme nearby
      # markers to zero for each contiguous region of extreme markers
      if (prune) {
        temp <- vector()
        for (i in 1:nrow(ret)) {
          if (! is.na(ret[i, ]$freq)) {
            temp <- c(temp, i)
          } else if (length(temp) > 0) {
            highest <- which(ret[temp, "freq"][[1]] == max(ret[temp, "freq"]))
            ret[temp[-highest], "freq"] <- NA
            temp <- vector()
          }
        }
      }
      print(ret)
      ret
    }
  ) %>% do.call(rbind, .)
}

load_groups <- function(csv) {
  read_csv(
    str_c("Data/Intermediate/Aligned_genes/selected_alignments/", csv),
    col_names = c("group", "chrom", "pos", "sleng", "salign", "%id")
  ) %>%
    mutate(pos_mb = pos / 1e6) %>%
    select(group, chrom, pos_mb) %>%
    cbind(base = 0)
}