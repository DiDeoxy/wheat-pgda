library(tidyverse)

calc_eh <- function (genotypes) {
  apply(genotypes, 1, function (snp) {
    2 * ((sum(snp == 0) / sum(snp == 0 | 2)) *
      (sum(snp == 2) / sum(snp == 0 | 2)))
  })
}

calc_max_genome_lengths <- function(wheat_data) {
  by(wheat_data$snp$pos_mb, wheat_data$snp$chrom, max) %>% 
    (
      function (max_chrom_lengths) {
        c(
          A = max(max_chrom_lengths[seq(1, 19, 3)]),
          B = max(max_chrom_lengths[seq(2, 20, 3)]),
          D = max(max_chrom_lengths[seq(3, 21, 3)])
        )
      }
    )
}

calc_max_homeolog_lengths <- function(wheat_data) {
  by(wheat_data$snp$pos, wheat_data$snp$chrom, max) %>%
    (
      function (max_chrom_lengths) {
        c(
          one = max(max_chrom_lengths[c(1, 2, 3)]),
          two = max(max_chrom_lengths[c(4, 5, 6)]),
          three = max(max_chrom_lengths[c(7, 8, 9)]),
          four = max(max_chrom_lengths[c(10, 11, 12)]),
          five = max(max_chrom_lengths[c(13, 14, 15)]),
          six = max(max_chrom_lengths[c(16, 17, 18)]),
          seven = max(max_chrom_lengths[c(19, 20, 21)])
        )
      }
    )
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

load_groups <- function(csv) {
  read_csv(
    str_c("Data/Intermediate/Aligned_genes/selected_alignments/", csv),
    col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
  ) %>%
    mutate(pos_mb = pos / 1e6) %>%
    select(id, chrom, pos_mb) %>%
    cbind(base = 0)
}

snp_densities <- function(chrom) {
  # the distances between markers on the chromosome
  gaps <- diff(chrom$pos_mb)
  # a vector for holding the average density of snps near snp i
  densities <- vector()
  for (i in 1:nrow(chrom)) {
    # a vector for holding the inices of the 10 snps upstream and 10 snps
    # downstream from snp i, as well as snp i
    nearest <- vector()
    if (i < 5 && length(gaps) >= 8) {
      # position 1 in gaps is the gap between snps 1 and 2
      # position 20 in gaps is the gap between snps 20 and 21
      nearest <- 1:(i + 4)
    } else if (i > length(gaps) - 5 && length(gaps) >= 8) {
      # position length(gaps) - 20 in gaps is the gap between snps
      # length(chrom) - 21  and length(chrom) - 20
      # position length(gaps) in gaps is the gap between snp length(chrom) - 1
      # and length(chrom)
      nearest <- (i - 4):length(gaps)
    } else if (length(gaps) < 8) {
      nearest <- 1:length(gaps)
    } else {
      # the positon i - 9 in gaps is the gap between snps i - 10 and i - 9
      # the psotion i + 9 in gaps is the gap between snps i + 9 and i + 10
      nearest <- (i - 4):(i + 4)
    }
    # calc the density of the regions aorund marker i and add it to the
    # densities vector
    densities <- c(densities, mean(gaps[nearest]))
  }
  densities
}

calc_group_extreme_freqs <- function(wheat_data, extremes) {
  by(wheat_data$snp, wheat_data$snp$chrom,
    # wheat_data$snp[which(wheat_data$snp$chrom == "5D"), ], wheat_data$snp$chrom[which(wheat_data$snp$chrom == "5D")],
    function(snp_data) {
      # calc the local densities of each marker
      snp_data <- snp_data %>% add_column(density = snp_densities(snp_data))
      # initialize an empty tibble with columns for the important data
      inter <- tibble(
        chrom = character(), group = character(), pos_mb = double(),
        freq = double(), num = integer(),
        extreme_D = double(),
        id = character()
      )
      # 
      for (group in names(extremes)) {
        for (i in 1:nrow(snp_data)) {
          nearby <- which(
            (snp_data$pos_mb >= 
              snp_data$pos_mb[i] - min((snp_data$density[i] * 5), 1)
            ) &
            (snp_data$pos_mb <= 
              snp_data$pos_mb[i] + min((snp_data$density[i] * 5), 1)
            )
          )
          freq <- sum(
            snp_data[[group]][nearby] > extremes[[group]]
          ) / length(nearby)
          inter <- inter %>%
            add_row(
              chrom = snp_data$chrom[i], group = group,
              pos_mb = snp_data$pos_mb[i], freq = freq, num = NA,
              extreme_D = ifelse(
                snp_data[[group]][i] > extremes[[group]], snp_data[[group]][i],
                NA
              ),
              id = snp_data$id[i]
            )
        }
      }
      temp <- vector()
      ret <- tibble(
        chrom = character(), group = character(), pos_mb = double(),
        freq = double(), num = integer(), extreme_Ds = character(),
        id = character()
      )
      ret <- ret %>% add_row(
        chrom = inter$chrom[i], group = unique(inter$group)
      )
      for (i in 1:nrow(inter)) {
        if (inter$freq[i] > 0) {
          temp <- c(temp, i)
        }
        if (((inter$freq[i] == 0 && length(temp) > 0) ||
             (i == nrow(inter) && length(temp) > 0)
            ) && mean(inter$freq[temp]) >= 0.2) {
              ret <- ret %>% add_row(
                chrom = inter$chrom[i], group = inter$group[i],
                pos_mb = mean(inter$pos_mb[temp]),
                freq = mean(inter$freq[temp]), num = length(temp), 
                extreme_Ds = paste(inter$extreme_D[temp], collapse= ' '),
                id = paste(inter$id[temp], collapse= ' ')
              )
              temp <- vector()
        }
      }
      print(ret, n = 50)
      ret
    }
  ) %>% do.call(rbind, .)
}