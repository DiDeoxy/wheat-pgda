library(tidyverse)
library(plyr)
library(GGally)
library(SNPRelate)
library(extrafont)

# load custom functions
source("src\\R_functions\\funcs_calc_stats.R")
source("src\\R_functions\\colour_sets.R")

# size of a megabase, used to divide the bp positions
mb <- 1000000

# load the data from the gds object
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")
snp_pos <- snp_pos / mb

# create a data frame of all the 
data_eh <- tibble(id = snp_id, chrom = snp_chrom, pos = snp_pos) %>%
           add_column(eh = eh(genotypes))

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(data_eh$pos, data_eh$chrom, max)
max_genome_lengths <- data.frame(A = max(chrom_lengths[seq(1, 19, 3)]),
                                 B = max(chrom_lengths[seq(2, 20, 3)]),
                                 D = max(chrom_lengths[seq(3, 21, 3)]))

# allows application of same colour to each set of chromosomes
colour_order <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3),
                  rep(7,3))

# identify those snps within the extended haplotypes
index_A1 <- which(snp_chrom == 1)
snps_A1 <- snp_id[index_A1][which(snp_pos[index_A1] > 70 & 
                                  snp_pos[index_A1] < 300)]
index_A2 <- which(snp_chrom == 4)
snps_A2 <- snp_id[index_A2][which(snp_pos[index_A2] > 210 &
                                  snp_pos[index_A2] < 470)]
index_A4 <- which(snp_chrom == 10)
snps_A4 <- snp_id[index_A4][which(snp_pos[index_A4] > 230 &
                                  snp_pos[index_A4] < 460)]
index_B5 <- which(snp_chrom == 14)
snps_B5 <- snp_id[index_B5][which(snp_pos[index_B5] > 110 &
                                  snp_pos[index_B5] < 210)]
index_A6 <- which(snp_chrom == 16)
snps_A6 <- snp_id[index_A6][which(snp_pos[index_A6] > 170 &
                                  snp_pos[index_A6] < 445)]
index_B6 <- which(snp_chrom == 17)
snps_B6 <- snp_id[index_B6][which(snp_pos[index_B6] > 250 &
                                  snp_pos[index_B6] < 380)]
index_A7 <- which(snp_chrom == 19)
snps_A7 <- snp_id[index_A7][which(snp_pos[index_A7] > 310 &
                                  snp_pos[index_A7] < 445)]
haplo_snps <- match(c(snps_A1, snps_A2, snps_A4, snps_B5, snps_A6, snps_B6,
                      snps_A7),
                    data$id)

# create a column in the data set that contains phi values for only those
# markers in the extended haplotypes, NA for all else
data_eh_haplo <- data_eh %>% mutate(haplo = eh)
data_eh_haplo$haplo[-haplo_snps] <- NA

# create plots of all markers' eh values on each chromosome
plots <- by(data_eh_haplo, data_eh_haplo$chrom, function (data_chrom) {
  chrom_num <- data_chrom$chrom[1]
  data_chrom %>%
    ggplot() +
      xlim(0, ifelse(chrom_num %in% seq(1, 19, 3), max_genome_lengths$A,
                      ifelse(chrom_num %in% seq(2, 30, 3),
                             max_genome_lengths$B, max_genome_lengths$D))) +
      geom_point(aes(pos, eh),
                 colour = colours_chroms[colour_order[chrom_num]],
                 size = 0.5) +
      geom_point(aes(pos, haplo), shape = 1, size = 0.75)
})

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "EH",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "EH of Markers"
)

# plot the matrix
png(str_c("Results\\loci\\EH\\full_EH.png"),
    family = "Times New Roman", width = 200, height = 287, pointsize = 5,
    units = "mm", res = 300)
plots_matrix
dev.off()