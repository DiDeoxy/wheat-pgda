library(tidyverse)
library(ggplot2)
library(GGally)
library(extrafont)
library(RColorBrewer)
library(pracma)

source("src/R_functions/funcs_gds_parse_create.R")

plot_gaps_nbs_ld <- function(lng, genome_ld, subset, plot_title) {
  # histograms and boxplots depicting the distribution of gaps on each genome
  gaps_log10 <- tibble(
    Genome = factor(
      c(rep("A", length(lng$A$gaps)),
        rep("B", length(lng$B$gaps)),
        rep("D", length(lng$D$gaps)),
        rep("All", length(c(lng$A$gaps, lng$B$gaps, lng$D$gaps)))
      ),
      levels = c("A", "B", "D", "All")
    ),
    gaps = c(
      lng$A$gaps,
      lng$B$gaps,
      lng$D$gaps,
      lng$A$gaps, lng$B$gaps, lng$D$gaps
    ) %>% log10()
  )

  # histograms and boxplots depicting the distribution of gaps on each genome
  nbs_ld_genome <- tibble(
    Genome = factor(
      c(rep("A", length(genome_ld$A$nbs)),
        rep("B", length(genome_ld$B$nbs)),
        rep("D", length(genome_ld$D$nbs)),
        rep("All", length(c(genome_ld$A$nbs, genome_ld$B$nbs, genome_ld$D$nbs)))
      ),
      levels = c("A", "B", "D", "All")
    ),
    ld = c(
      genome_ld$A$nbs,
      genome_ld$B$nbs,
      genome_ld$D$nbs,
      genome_ld$A$nbs, genome_ld$B$nbs, genome_ld$D$nbs
    )
  )

  plots <- list()
  plots[[1]] <- gaps_log10 %>%
    ggplot() +
    geom_freqpoly(aes(gaps, colour = Genome), size = 0.3) +
    scale_color_manual(values = brewer.pal(4, "Dark2")) +
    xlim(min(gaps_log10$gaps), max(gaps_log10$gaps)) +
    ylim(0, 500)
  plots[[2]] <- nbs_ld_genome %>%
    ggplot() +
    geom_freqpoly(aes(ld, colour = Genome), size = 0.3) +
    scale_color_manual(values = brewer.pal(4, "Dark2")) +
    xlim(0, 1.0001) +
    ylim(0, 500)

  # turn plot list into ggmatrix
  plots_matrix <- ggmatrix(
    plots, nrow = 1, ncol = 2, yAxisLabels = "Num Markers",
    xAxisLabels = c(
      "Log 10 of Mb Gap Distances",
      "Absolute Composite LD of Neighbouring Markers"
    ),
    title = plot_title,
    legend = c(1, 2)
  )

  # plot the matrix
  png(str_c("Results/gaps/gaps_nbs_", subset, ".png"),
    family = "Times New Roman", width = 100, height = 62, pointsize = 10,
    units = "mm", res = 300)
  print(plots_matrix + 
    theme(legend.position = "bottom",
      text = element_text(
        size = 8, lineheight = 0.1)))
  dev.off()
}

calc_ld_stats <- function (subset, snp_data) {
  wheat_internal <- snpgdsOpen(
    str_c("Data/Intermediate/GDS/", subset, ".gds"))
  # Calcualte ld between all snps on each chromosome
  ld_stats <- by(snp_data, snp_data$chrom, function (chrom) {
    ld_mat <- snpgdsLDMat(wheat_internal, method = "composite",
      slide = -1, snp.id = chrom$id)$LD %>%
      abs() 
    return(
      list(
        pw = ld_mat %>% as.dist() %>% as.vector(),
        nbs = ld_mat %>% Diag(., 1)
      )
    )
  })
  snpgdsClose(wheat_internal)

  genome_ld <- list(A = list(), B = list(), D = list())
  for (i in 1:length(ld_stats)) {
    if (i %in% seq(1, 19, 3)) {
      genome_ld$A$pw <- c(genome_ld$A$pw, ld_stats[[i]]$pw)
      genome_ld$A$nbs <- c(genome_ld$A$nbs, ld_stats[[i]]$nbs)
    } else if (i %in% seq(2, 20, 3)) {
      genome_ld$B$pw <- c(genome_ld$B$pw, ld_stats[[i]]$pw)
      genome_ld$B$nbs <- c(genome_ld$B$nbs, ld_stats[[i]]$nbs)
    } else {
      genome_ld$D$pw <- c(genome_ld$D$pw, ld_stats[[i]]$pw)
      genome_ld$D$nbs <- c(genome_ld$D$nbs, ld_stats[[i]]$nbs)
    }
  }

  genome_ld
}

calc_maf_mr <- function (wheat_data) {
  chrom_geno_sums <- by(wheat_data$genotypes, wheat_data$snp$chrom,
    function (chrom_genos) {
      apply(chrom_genos, 1, function (snp) {
        A <- sum(snp == 0)
        B <- sum(snp == 2)
        missing <- sum(snp == 3)
        return(tibble(maf = min(c(A, B)) / sum(A, B), mr = missing / sum(A, B, missing)))
      }) %>% bind_rows()
    }
  )
  list(
    A = chrom_geno_sums[seq(1, 19, 3)] %>% bind_rows(),
    B = chrom_geno_sums[seq(2, 20, 3)] %>% bind_rows(),
    D = chrom_geno_sums[seq(3, 21, 3)] %>% bind_rows()
  )
}

calc_lng <- function(snp_data) {
  # number of snps and mean distances between the genome
  lng <- by(snp_data, snp_data$chrom, 
    function (chrom) {
      list(
        leng = max(chrom$pos_mb),
        num = length(chrom$pos_mb),
        gaps = diff(chrom$pos_mb)
      )
    }
  )
  ret <- list(
    A = list(leng = vector(), num = vector(), gaps = vector()),
    B = list(leng = vector(), num = vector(), gaps = vector()),
    D = list(leng = vector(), num = vector(), gaps = vector())
  )
  for (i in 1:length(lng)) {
    if (i %in% seq(1, 19, 3)) {
      ret$A$leng <- c(ret$A$leng, lng[[i]]$leng)
      ret$A$num <- c(ret$A$num, lng[[i]]$num)
      ret$A$gaps <- c(ret$A$gaps, lng[[i]]$gaps)
    } else if (i %in% seq(2, 20, 3)) {
      ret$B$leng <- c(ret$B$leng, lng[[i]]$leng)
      ret$B$num <- c(ret$B$num, lng[[i]]$num)
      ret$B$gaps <- c(ret$B$gaps, lng[[i]]$gaps)
    } else {
      ret$D$leng <- c(ret$D$leng, lng[[i]]$leng)
      ret$D$num <- c(ret$D$num, lng[[i]]$num)
      ret$D$gaps <- c(ret$D$gaps, lng[[i]]$gaps)
    }
  }
  ret
}

calc_plot_map_stats <- function (subset, plot_title_1, plot_title_2) {
  wheat_data <- parse_gds(subset)

  # calc ld stats
  genome_ld <- calc_ld_stats(subset, wheat_data$snp)

  # find the most distant snp on each chroms, the number of snps on each,
  # and the sizes of the gaps between snps
  lng <- calc_lng(wheat_data$snp)

  # plot the ld
  plot_gaps_nbs_ld(lng, genome_ld, subset, plot_title_2)

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
  write_csv(map_stats, str_c("Results/gaps/map_stats_", subset, ".csv"))
}