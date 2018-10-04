library(tidyverse)
library(plyr)
library(GGally)
library(SNPRelate)
library(extrafont)

# load custom functions
source("src\\R_functions\\funcs_calc_stats.R")
source("src\\R_functions\\colour_sets.R")

plot_eh <- function (wheat, subset) {

  # create a snp_data frame of all the
  wheat$snp_data <- wheat$snp_data %>%
    add_column(eh = calc_eh(wheat$genotypes))

  # find the max position of any marker on each genome for xlims
  chrom_lengths <- by(wheat$snp_data$pos, wheat$snp_data$chrom, max)
  max_genome_lengths <- c(max(chrom_lengths[seq(1, 19, 3)]), # A genome
                          max(chrom_lengths[seq(2, 20, 3)]), # B genome
                          max(chrom_lengths[seq(3, 21, 3)])) # D genome

  # allows application of same colour to each set of chromosomes
  colour_order <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), 
    rep(6, 3), rep(7, 3))

  # identify those snps within the extended haplotypes
  index_A1 <- which(wheat$snp_data$chrom == 1)
  snps_A1 <- wheat$snp_data$id[index_A1][which(wheat$snp_data$pos[index_A1] > 70 &
                                    wheat$snp_data$pos[index_A1] < 300)]
  index_A2 <- which(wheat$snp_data$chrom == 4)
  snps_A2 <- wheat$snp_data$id[index_A2][which(wheat$snp_data$pos[index_A2] > 210 &
                                    wheat$snp_data$pos[index_A2] < 470)]
  index_A4 <- which(wheat$snp_data$chrom == 10)
  snps_A4 <- wheat$snp_data$id[index_A4][which(wheat$snp_data$pos[index_A4] > 230 &
                                    wheat$snp_data$pos[index_A4] < 460)]
  index_B5 <- which(wheat$snp_data$chrom == 14)
  snps_B5 <- wheat$snp_data$id[index_B5][which(wheat$snp_data$pos[index_B5] > 110 &
                                    wheat$snp_data$pos[index_B5] < 210)]
  index_A6 <- which(wheat$snp_data$chrom == 16)
  snps_A6 <- wheat$snp_data$id[index_A6][which(wheat$snp_data$pos[index_A6] > 170 &
                                    wheat$snp_data$pos[index_A6] < 445)]
  index_B6 <- which(wheat$snp_data$chrom == 17)
  snps_B6 <- wheat$snp_data$id[index_B6][which(wheat$snp_data$pos[index_B6] > 250 &
                                    wheat$snp_data$pos[index_B6] < 380)]
  index_A7 <- which(wheat$snp_data$chrom == 19)
  snps_A7 <- wheat$snp_data$id[index_A7][which(wheat$snp_data$pos[index_A7] > 310 &
                                    wheat$snp_data$pos[index_A7] < 445)]
  haplo_snps <- match(c(snps_A1, snps_A2, snps_A4, snps_B5, snps_A6, snps_B6,
                        snps_A7),
                      wheat$snp_data$id)

  # create a column in the snp_data set that contains phi values for only
  # those markers in the extended haplotypes, NA for all else
  wheat$snp_data <- wheat$snp_data %>% mutate(haplo = eh)
  wheat$snp_data[-haplo_snps] <- NA

  # create plots of all markers' eh values on each chromosome
  plots <- by(wheat$snp_data, wheat$snp_data$chrom, 
    function (snp_data_chrom) {
      chrom_num <- snp_data_chrom$chrom[1]
      snp_data_chrom %>%
        ggplot() +
          xlim(0, max_genome_lengths[ifelse(chrom_num %% 3, 
            chrom_num %% 3, 3)]) +
          geom_point(aes(pos, eh),
                    colour = colours_chroms[colour_order[chrom_num]],
                    size = 0.5) +
          geom_point(aes(pos, haplo), shape = 1, size = 0.75)
    })

  # turn plot list into ggmatrix
  plots_matrix <- ggmatrix(
    plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "EH",
    xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
    title = str_c("EH of Markers ", subset)
  )

  # plot the matrix
  png(str_c("Results\\loci\\EH\\", subset, "_EH.png"),
      family = "Times New Roman", width = 200, height = 287, pointsize = 5,
      units = "mm", res = 300)
  print(plots_matrix)
  dev.off()
}