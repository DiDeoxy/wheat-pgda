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

load_groups <- function(csv, base = 0) {
  read_csv(
    str_c("Data/Intermediate/Aligned_genes/selected_alignments/", csv),
    col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
  ) %>%
    mutate(pos_mb = pos / 1e6) %>%
    select(id, chrom, pos_mb) %>%
    cbind(base = base)
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
    function(snp_data) {
      # initialize an empty tibble with the needed columns
      inter <- tibble(
        chrom = character(), pos_mb = double(),
        num_nearby_extreme = double(), D = double(), id = character()
      )
      groups <- list("chrs_csws" = inter, "chrs_chrw" = inter, "csws_chrw" = inter)
      # the number of markers upstream and downstream of the one considered for
      # the identification of the local frequency of extreme markers
      # num <- 10
      # the max distance in Mb from the considered marker that the upstream and
      # downstream markers of num can be
      dist <- 4
      # for each comparison
      for (group in names(extremes)) {
        # for each snp in the comparison
        for (i in 1:nrow(snp_data)) {
          # for holding those markers considered nearby (within num & dist)
          nearby <- i
          # find upstream nearby markers
          if (i > 1) {
            for (j in (i - 1):1) {
              if (
                j >= 1
                && snp_data$pos_mb[i] - snp_data$pos_mb[j] <= dist
                # && i - j <= num
              ) {
                nearby <- c(j, nearby)
              } else {
                break
              }
            }
          }
          # find downstream nearby markers
          if (i < nrow(snp_data)) {
            for (j in (i + 1):nrow(snp_data)) {
              if (
                j <= nrow(snp_data)
                && snp_data$pos_mb[j] - snp_data$pos_mb[i] <= dist
                # && j - i <= num
              ) {
                nearby <- c(nearby, j)
              } else {
                break
              }
            }
          }
          # store the useful data for each marker in each group including its
          # frequency of nearby extreme markers, its Jost's D values, and its id
          groups[[group]] <- groups[[group]] %>%
            add_row(
              chrom = snp_data$chrom[i], pos_mb = snp_data$pos_mb[i],
              num_nearby_extreme = sum(
                snp_data[[group]][nearby] > extremes[[group]]
              ),
              D = snp_data[[group]][i], id = snp_data$id[i]
            )
        }
      }

      linked <- vector()
      ret <- tibble(
        chrom = character(), group = character(), mean_pos_mb = double(),
        num_linked = integer(), freq_extreme = double(), mean_D = double(),
        extreme = character(), pos_mb = character(), Ds = character(),
        ids = character()
      )
      ret <- ret %>% add_row(
        chrom = groups[[group]]$chrom[i], group = names(groups)
      )
      for (group in names(groups)) {
        for (row in 1:nrow(groups[[group]])) {
          if (groups[[group]]$num_nearby_extreme[row]) {
            # print(groups[[group]]$num_nearby_extreme[row])
            linked <- c(linked, row)
          }
          # needs to be a separate else because the last row can be in linked
          if (
            (
              ! groups[[group]]$num_nearby_extreme[row] ||
              row == nrow(groups[[group]])
            )
            && length(linked)
          ) {
            linked_extreme <- which(groups[[group]]$D[linked] > extremes[group])
            linked_pruned <- linked[
              linked_extreme[1]:linked_extreme[length(linked_extreme)]
            ]
            pos_pruned <- groups[[group]]$pos_mb[linked_pruned]
            if (
              # when there are markers within dist to a single extreme marker,
              # skip it
              (length(linked) > 1 && length(linked_extreme) == 1) ||
              # when there are maore than 1 extreme marker in a region less than
              # 1.5 Mb with distal markers that are not extreme and more than
              # 75% of them are extreme skip it (looks like a bunch of markers
              # acting like a single marker)
              (
                length(linked_pruned) > 1
                && pos_pruned[length(pos_pruned)] - pos_pruned[1] <= 1.5
                && length(linked_pruned) < length(linked)
                && length(linked_extreme) / length(linked_pruned) > 0.75
              )
            ) {
              linked <- vector()
              next
            }
            if (
              mean(groups[[group]]$D[linked_pruned]) >= 0
            ) {
              ret <- ret %>% add_row(
                chrom = groups[[group]]$chrom[i], group = group,
                mean_pos_mb = mean(groups[[group]]$pos_mb[linked_pruned]),
                num_linked = length(linked),
                freq_extreme =
                  sum(
                    groups[[group]]$D[linked_pruned] >= extremes[group]
                  ) / length(linked_pruned),
                mean_D = mean(groups[[group]]$D[linked_pruned]),
                extreme = paste(
                  groups[[group]]$D[linked_pruned] > extremes[group],
                  collapse = ' '
                ),
                pos_mb = paste(
                  groups[[group]]$pos_mb[linked_pruned], collapse = ' '
                ),
                Ds = paste(groups[[group]]$D[linked_pruned], collapse = ' '),
                ids = paste(groups[[group]]$id[linked_pruned], collapse = ' ')
              )
            }
            linked <- vector()
          }
        }
      }
      print(ret, n = 50)
      ret
    }
  ) %>% do.call(rbind, .)
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits) * 10^digits

  (df)
}