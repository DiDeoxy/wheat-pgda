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

plot_neighbour_ld <- function(subset, snp_data, plot_title) {
  wheat_internal <- snpgdsOpen(
    str_c("Data/Intermediate/GDS/", subset, ".gds"))
  # Calcualte ld between all snps on each chromosome
  neighbour_ld_chrom <- by(snp_data, snp_data$chrom, function (chrom) {
    neighbour_ld <- snpgdsLDMat(wheat_internal, method = "composite",
      slide = -1, snp.id = chrom$id)$LD %>%
      abs() %>%
      Diag(., 1)
  })
  snpgdsClose(wheat_internal)

  A <- unlist(neighbour_ld_chrom[seq(1, 19, 3)])
  B <- unlist(neighbour_ld_chrom[seq(2, 20, 3)])
  D <- unlist(neighbour_ld_chrom[seq(3, 21, 3)])

  # histograms and boxplots depicting the distribution of gaps on each genome
  neighbour_ld_genome <- tibble(
    Genome = factor(
      c(rep("A", length(A)),
        rep("B", length(B)),
        rep("D", length(D)),
        rep("All", length(c(A, B, D)))
      ),
      levels = c("A", "B", "D", "All")
    ),
    ld = c(
      A, B, D,
      A, B, D
    )
  )

  plots <- list()
  plots[[2]] <- ggplot(neighbour_ld_genome, aes(ld, colour = Genome)) +
    geom_freqpoly() +
    scale_color_manual(values = brewer.pal(4, "Dark2")) +
    xlim(0, 1)
  neighbour_ld_genome$Genome <- factor(neighbour_ld_genome$Genome,
    rev(levels(neighbour_ld_genome$Genome)))
  plots[[1]] <- ggplot(neighbour_ld_genome, aes(Genome, ld, colour = Genome)) +
    geom_boxplot() +
    ylim(0, 1) +
    coord_flip() +
    scale_color_manual(values = rev(brewer.pal(4, "Dark2")))

  # turn plot list into ggmatrix
  plots_matrix <- ggmatrix(
    plots, nrow = 2, ncol = 1, xlab = "LD between Neighbours",
    yAxisLabels = c("a) Boxplots", "b) Frequency Plots"),
    title = plot_title,
    legend = c(2, 1)
  )

  # plot the matrix
  png(str_c("Results/gaps/neighbour_ld_", subset, ".png"),
    family = "Times New Roman", width = 100, height = 143, pointsize = 10,
    units = "mm", res = 300)
  print(plots_matrix + theme(legend.position = "bottom"))
  dev.off()
}

plot_gaps <- function(gaps, subset, plot_title) {
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

  plots <- list()
  plots[[2]] <- ggplot(gaps_log10, aes(gaps, colour = Genome)) +
    geom_freqpoly() +
    scale_color_manual(values = brewer.pal(4, "Dark2")) +
    xlim(min(gaps_log10$gaps), max(gaps_log10$gaps))
  gaps_log10$Genome <- factor(gaps_log10$Genome,
    rev(levels(gaps_log10$Genome)))
  plots[[1]] <- ggplot(gaps_log10, aes(Genome, gaps, colour = Genome)) +
    geom_boxplot() +
    ylim(min(gaps_log10$gaps), max(gaps_log10$gaps)) +
    coord_flip() +
    scale_color_manual(values = rev(brewer.pal(4, "Dark2")))

  # turn plot list into ggmatrix
  plots_matrix <- ggmatrix(
    plots, nrow = 2, ncol = 1, xlab = "Log10 Transformed Gap Distance",
    yAxisLabels = c("a) Boxplots", "b) Frequency Plots"),
    title = plot_title,
    legend = c(2, 1)
  )

  # plot the matrix
  png(str_c("Results/gaps/gaps_dist_", subset, ".png"),
    family = "Times New Roman", width = 100, height = 143, pointsize = 10,
    units = "mm", res = 300)
  print(plots_matrix + theme(legend.position = "bottom"))
  dev.off()
}

calc_mean_pairwise_ld <- function (subset, snp_data) {
  wheat_internal <- snpgdsOpen(
    str_c("Data/Intermediate/GDS/", subset, ".gds"))
  # Calcualte ld between all snps on each chromosome
  mean_ld <- by(snp_data, snp_data$chrom, function (chrom) {
    snpgdsLDMat(wheat_internal, method = "composite",
      slide = -1, snp.id = chrom$id)$LD %>%
      as.dist %>% abs() %>% mean(na.rm = TRUE) 
  })
  snpgdsClose(wheat_internal)
  mean_ld
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

  # plot the gaps
  plot_gaps(gaps, subset, plot_title_1)

  # plot the ld
  plot_neighbour_ld(subset, wheat_data$snp, plot_title_2)

  # find the min length of the top percentile of gaps
  top_percentile <- quantile(unlist(gaps), prob = 0.99, na.rm = T)

  mean_pairwise_ld <- calc_mean_pairwise_ld(subset, wheat_data$snp)

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
    mean(unlist(mean_pairwise_ld[seq(1, 19, 3)]), na.rm = TRUE), "\n",
    "Genome B: average pairwise ld ",
    mean(unlist(mean_pairwise_ld[seq(2, 20, 3)]), na.rm = TRUE), "\n",
    "Genome D: average pairwise ld ",
    mean(unlist(mean_pairwise_ld[seq(3, 21, 3)]), na.rm = TRUE), "\n",
    "Overall average pairwise ld ",
    mean(unlist(mean_pairwise_ld), na.rm = TRUE), "\n"
  )
  out <- file(str_c("Results/gaps/gaps_dist_", subset, ".txt"))
  writeLines(report, out)
  close(out)
}