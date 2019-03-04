#' Calculate and plot map stats
#'
#' Calcluates various map statistics creating plots and tables
#'
#' @importFrom dplyr tibble
#' @importFrom readr write_csv
#' @importFrom stringr str_c
#'
#' @param gds the file path to the gds object
#' @param plot_title the title of the combined plots
#' @param y_lim the limit for the y axis of the plots
#' @param out_name the base name to use for the csv and png outputs
#'
#' @return Outputs plots and tables of the data
#'
#' @export

calc_plot_map_stats <- function (gds, plot_title, y_lim, out_name) {
  wheat_data <- snpgds_parse(gds)

  # calc ld stats
  genome_ld <- calc_ld_stats(gds, wheat_data$snp)

  # find the most distant snp on each chroms, the number of snps on each,
  # and the sizes of the gaps between snps
  lng <- calc_lng(wheat_data$snp)

  # plot the ld
  plot_gaps_nbs_ld(lng, genome_ld, gds, plot_title, y_lim, out_name)

  # find the min length of the top percentile of gaps
  top_percentile <- quantile(
    c(lng$A$gaps, lng$B$gaps, lng$D$gaps),
    prob = 0.99, na.rm = T
  )

  # find the maf and mr
  maf_mr <- calc_maf_mr(wheat_data)

  map_stats <- tibble(
    "Genome" = c("A", "B", "D", "All"),
    "MAF" = c(
      mean(maf_mr$A$maf),
      mean(maf_mr$B$maf),
      mean(maf_mr$D$maf),
      mean(c(maf_mr$A$maf, maf_mr$B$maf, maf_mr$D$maf))
    ),
    "MR" = c(
      mean(maf_mr$A$mr),
      mean(maf_mr$B$mr),
      mean(maf_mr$D$mr),
      mean(c(maf_mr$A$mr, maf_mr$B$mr, maf_mr$D$mr))
    ),
    "Covered Bases (Gb)" = c(
      sum(lng$A$leng) / 1000,
      sum(lng$B$leng) / 1000,
      sum(lng$D$leng) / 1000,
      sum(lng$A$leng, lng$B$leng, lng$D$leng) / 1000
    ),
    "Num. SNPs" = c(
      sum(lng$A$num),
      sum(lng$B$num),
      sum(lng$D$num),
      sum(lng$A$num, lng$B$num, lng$D$num)
    ),
    "Mean Gap Size (Mb)" = c(
      mean(lng$A$gaps),
      mean(lng$B$gaps),
      mean(lng$D$gaps),
      mean(c(lng$A$gaps, lng$B$gaps, lng$D$gaps))
    ),
    "Num. Top 1% Gaps" = c(
      sum(lng$A$gaps >= top_percentile),
      sum(lng$B$gaps >= top_percentile),
      sum(lng$D$gaps >= top_percentile),
      sum(
        c(
          lng$A$gaps >= top_percentile,
          lng$B$gaps >= top_percentile,
          lng$D$gaps >= top_percentile
        )
      )
    ),
    "Min Length Top 1% Gap (Mb)" = c(
      min(lng$A$gaps[which(lng$A$gaps >= top_percentile)]),
      min(lng$B$gaps[which(lng$B$gaps >= top_percentile)]),
      min(lng$D$gaps[which(lng$D$gaps >= top_percentile)]),
      top_percentile
    ),
    "Max Length Gap (Mb)" = c(
      max(lng$A$gaps),
      max(lng$B$gaps),
      max(lng$D$gaps),
      max(c(lng$A$gaps, lng$B$gaps, lng$D$gaps))
    ),
    "Mean Pairwise LD" = c(
      mean(genome_ld$A$pw, na.rm = TRUE),
      mean(genome_ld$B$pw, na.rm = TRUE),
      mean(genome_ld$D$pw, na.rm = TRUE),
      mean(c(genome_ld$A$pw, genome_ld$B$pw, genome_ld$D$pw), na.rm = TRUE)
    ),
    "Mean Neighbour LD" = c(
      mean(genome_ld$A$nbs, na.rm = TRUE),
      mean(genome_ld$B$nbs, na.rm = TRUE),
      mean(genome_ld$D$nbs, na.rm = TRUE),
      mean(c(genome_ld$A$nbs, genome_ld$B$nbs, genome_ld$D$nbs), na.rm = TRUE)
    )
  )
  write_csv(
    map_stats, file.path(
      map_stats_and_plots, str_c(out_name, ".csv")
    )
  )
}