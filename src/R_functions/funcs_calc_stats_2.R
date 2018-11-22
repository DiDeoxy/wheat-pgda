library(tidyverse)
library(ggplot2)
library(GGally)
library(extrafont)
library(RColorBrewer)
# install.packages("pracma")
library(pracma)

source("src/R_functions/funcs_gds_parse_create.R")

calc_eh <- function (genotypes) {
  apply(genotypes, 1, function (snp) {
    2 * ((sum(snp == 0) / sum(snp == 0 | 2)) *
      (sum(snp == 2) / sum(snp == 0 | 2)))
  })
}

calc_maf_mr <- function (genotypes) {
  # count how often for each snp each homozygote occurs
  counts <- apply(genotypes, 1, function (snp) {
    count <- as.data.frame(table(snp))
    if (nrow(count) == 3) {
      return(count[, 2])
    } else {
      return(c(count[, 2], 0))
    }
  })
  # find the minor allele frequency and missing rate of each snp
  as.data.frame(t(apply(counts, 2, function (snp) {
    maf <- min(snp[1:2]) / sum(snp[1:2])
    mr <- snp[3] / sum(snp)
    return(c(maf = maf, mr = mr))
  })))
}

plot_gaps_nbs_ld <- function(gaps, genome_ld, subset, plot_title) {
  # histograms and boxplots depicting the distribution of gaps on each genome
  gaps_log10 <- tibble(
    Genome = factor(
      c(rep("A", length(unlist(gaps$A))),
        rep("B", length(unlist(gaps$B))),
        rep("D", length(unlist(gaps$D))),
        rep("All", length(unlist(gaps)))
      ),
      levels = c("A", "B", "D", "All")
    ),
    gaps = log10(
      c(unlist(gaps$A), unlist(gaps$B), unlist(gaps$D), unlist(gaps))
    )
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
    xlim(0, 1) +
    ylim(0, 500)

  # turn plot list into ggmatrix
  plots_matrix <- ggmatrix(
    plots, nrow = 1, ncol = 2, yAxisLabels = "Num Markers",
    xAxisLabels = c("Gap Distances", "Neigbouring LD"),
    title = plot_title,
    legend = c(1, 2)
  )

  # plot the matrix
  png(str_c("Results/gaps/gaps_nbs_ld_", subset, ".png"),
    family = "Times New Roman", width = 100, height = 62, pointsize = 10,
    units = "mm", res = 300)
  print(plots_matrix + 
    theme(legend.position = "bottom",
      text = element_text(
        size = 8, lineheight = 0.1)))
  dev.off()
}
theme()
# subset <- "phys_subset_sample_ld_pruned"
# wheat_data <- parse_gds(subset)
# snp_data <- wheat_data$snp
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

calc_length_num_gaps <- function(snp_data) {
  ## number of snps and mean distances between the genome
  leng <- list(A = vector(), B = vector(), D = vector())
  num <- list(A = vector(), B = vector(), D = vector())
  gaps <- list(A = list(), B = list(), D = list())
  blah <- by(snp_data, snp_data$chrom, function (chrom) {
    if (chrom$chrom[1] %in% seq(1, 19, 3)) {
      leng$A <<- c(leng$A, max(chrom$pos_mb))
      num$A <<- c(num$A, length(chrom$pos_mb)) # number of snps on group
      gaps$A[[chrom$chrom[1]]] <<- diff(chrom$pos_mb)
    } else if (chrom$chrom[1] %in% seq(2, 20, 3)) {
      leng$B <<- c(leng$B, max(chrom$pos_mb))
      num$B <<- c(num$B, length(chrom$pos_mb))
      gaps$B[[chrom$chrom[1]]] <<- diff(chrom$pos_mb)
    } else {
      leng$D <<- c(leng$D, max(chrom$pos_mb))
      num$D <<- c(num$D, length(chrom$pos_mb))
      gaps$D[[chrom$chrom[1]]] <<- diff(chrom$pos_mb)
    }
  })
  return(list(leng = leng, num = num, gaps = gaps))
}

calc_plot_map_stats <- function (subset, plot_title_1, plot_title_2) {
  wheat_data <- parse_gds(subset)
  # find the maf and mr
  maf_mr <- calc_maf_mr(wheat_data$genotypes)

  # fidn the most distant snp on each chroms, the number of snps on each,
  # and the sizes of the gaps between snps
  length_num_gaps <- calc_length_num_gaps(wheat_data$snp)
  leng <- length_num_gaps$leng
  num <- length_num_gaps$num
  gaps <- length_num_gaps$gaps

  # calc ld stats
  genome_ld <- calc_ld_stats(subset, wheat_data$snp)

  # plot the ld
  plot_gaps_nbs_ld(gaps, genome_ld, subset, plot_title_2)

  # find the min length of the top percentile of gaps
  top_percentile <- quantile(unlist(gaps), prob = 0.99, na.rm = T)

  # find which chrom of each genome has the longest gap
  chr_A <- str_c(lapply(gaps$A, max) %>% which.max() %/% 3 + 1, "A")
  chr_B <- str_c(lapply(gaps$B, max) %>% which.max() %/% 3 + 1, "B")
  chr_D <- str_c(lapply(gaps$D, max) %>% which.max() %/% 3 + 1, "D")

  # print out all the data nicely
  report <- str_c(
    "A maf of ", mean(maf_mr$maf), "\n",
    "A missing rate of ", mean(maf_mr$mr), "\n",
    "############################################\n",
    "Genome A: num SNPs ", sum(num$A), " covering ", sum(leng$A), " Mb.\n",
    "Genome B: num SNPs ", sum(num$B), " covering ", sum(leng$B), " Mb.\n",
    "Genome D: num SNPs ", sum(num$D), " covering ", sum(leng$D), " Mb.\n",
    "Overall: num SNPs ", sum(unlist(num)), " covering ", sum(unlist(leng)),
    "\n",
    "############################################\n",
    "Genome A: average gap size ", mean(unlist(gaps$A)), " Mb.\n",
    "Genome B: average gap size ", mean(unlist(gaps$B)), " Mb.\n",
    "Genome D: average gap size ", mean(unlist(gaps$D)), " Mb.\n",
    "Overall: average gap size ", mean(unlist(gaps)), " Mb.\n",
    "############################################\n",
    "Genome A: num top 1% gaps ", sum(unlist(gaps$A) >= top_percentile), "\n",
    "Genome B: num top 1% gaps ", sum(unlist(gaps$B) >= top_percentile), "\n",
    "Genome D: num top 1% gaps ", sum(unlist(gaps$D) >= top_percentile), "\n",
    "############################################\n",
    "Min length top 1% gap: ", top_percentile, "\n",
    "Genome A: longest gap ", max(unlist(gaps$A)), " on chr ", chr_A, ".\n",
    "Genome B: longest gap ", max(unlist(gaps$B)), " on chr ", chr_B, ".\n",
    "Genome D: longest gap ", max(unlist(gaps$D)), " on chr ", chr_D, ".\n",
    "############################################\n",
    "Genome A: average pairwise ld ",
    mean(genome_ld$A$pw, na.rm = TRUE), "\n",
    "Genome B: average pairwise ld ",
    mean(genome_ld$B$pw, na.rm = TRUE), "\n",
    "Genome D: average pairwise ld ",
    mean(genome_ld$D$pw, na.rm = TRUE), "\n",
    "Overall average pairwise ld ",
    mean(c(genome_ld$A$pw, genome_ld$B$pw, genome_ld$D$pw), na.rm = TRUE), "\n",
    "############################################\n",
    "Genome A: average neighbouring ld ",
    mean(genome_ld$A$nbs), "\n",
    "Genome B: average neighbouring ld ",
    mean(genome_ld$B$nbs), "\n",
    "Genome D: average neighbouring ld ",
    mean(genome_ld$D$nbs), "\n",
    "Overall average neighbouring ld ",
    mean(c(genome_ld$A$nbs, genome_ld$B$nbs, genome_ld$D$nbs)), "\n"
  )
  out <- file(str_c("Results/gaps/gaps_dist_", subset, ".txt"))
  writeLines(report, out)
  close(out)
}