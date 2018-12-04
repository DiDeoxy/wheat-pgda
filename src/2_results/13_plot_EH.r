library(tidyverse)
library(plyr)
library(GGally)
library(SNPRelate)
library(extrafont)

# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/funcs_calc_stats.R")
source("src/R_functions/colour_sets.R")

# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

# create a snp_data frame of all the
wheat_data$snp <- wheat_data$snp %>%
  add_column(eh = calc_eh(wheat_data$genotypes))

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(wheat_data$snp$pos_mb, wheat_data$snp$chrom, max)
max_genome_lengths <- c(
  max(chrom_lengths[seq(1, 19, 3)]), # A genome
  max(chrom_lengths[seq(2, 20, 3)]), # B genome
  max(chrom_lengths[seq(3, 21, 3)])  # D genome
)

# identify those snps within the extended haplotypes
snp_index_1A <- which(wheat_data$snp$chrom == 1)
haplo_id_1A <- wheat_data$snp$id[snp_index_1A][
  which(wheat_data$snp$pos_mb[snp_index_1A] > 70 &
    wheat_data$snp$pos_mb[snp_index_1A] < 300)
]
snp_index_2A <- which(wheat_data$snp$chrom == 4)
haplo_id_2A <- wheat_data$snp$id[snp_index_2A][
  which(wheat_data$snp$pos_mb[snp_index_2A] > 210 &
    wheat_data$snp$pos_mb[snp_index_2A] < 470)
]
snp_index_4A <- which(wheat_data$snp$chrom == 10)
haplo_id_4A <- wheat_data$snp$id[snp_index_4A][
  which(wheat_data$snp$pos_mb[snp_index_4A] > 230 &
    wheat_data$snp$pos_mb[snp_index_4A] < 460)
]
snp_index_5B <- which(wheat_data$snp$chrom == 14)
haplo_id_5B <- wheat_data$snp$id[snp_index_5B][
  which(wheat_data$snp$pos_mb[snp_index_5B] > 110 &
    wheat_data$snp$pos_mb[snp_index_5B] < 210)
]
snp_index_6A <- which(wheat_data$snp$chrom == 16)
haplo_id_6A <- wheat_data$snp$id[snp_index_6A][
  which(wheat_data$snp$pos_mb[snp_index_6A] > 170 &
    wheat_data$snp$pos_mb[snp_index_6A] < 445)
]
snp_index_6B <- which(wheat_data$snp$chrom == 17)
haplo_id_6B <- wheat_data$snp$id[snp_index_6B][
  which(wheat_data$snp$pos_mb[snp_index_6B] > 250 &
    wheat_data$snp$pos_mb[snp_index_6B] < 380)
]
snp_index_7A <- which(wheat_data$snp$chrom == 19)
haplo_id_7A <- wheat_data$snp$id[snp_index_7A][
  which(wheat_data$snp$pos_mb[snp_index_7A] > 310 &
    wheat_data$snp$pos_mb[snp_index_7A] < 445)
]
haplo_ids <- c(
  haplo_id_1A, haplo_id_2A, haplo_id_4A, haplo_id_5B, haplo_id_6A, haplo_id_6B,
  haplo_id_7A
)
haplo_index_snps <- match(haplo_ids, wheat_data$snp$id)

# create a column in the snp_data set that contains phi values for only
# those markers in the extended haplotypes, NA for all else
wheat_data$snp <- wheat_data$snp %>% mutate(haplo = eh)
wheat_data$snp$haplo[-haplo_index_snps] <- NA

# allows application of same colour to each set of chromosomes
colour_order <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3),
rep(6, 3), rep(7, 3))

# create plots of all markers' eh values on each chromosome
plots <- by(wheat_data$snp, wheat_data$snp$chrom,
  function (data_chrom) {
    chrom_num <- data_chrom$chrom[1]
    data_chrom %>%
      ggplot() +
      xlim(
        0,
        max_genome_lengths[ifelse(chrom_num %% 3, chrom_num %% 3, 3)]
      ) +
      geom_point(
        aes(pos_mb, eh), colour = colours_chroms[colour_order[chrom_num]],
        size = 0.5
      ) +
      geom_point(aes(pos_mb, haplo), shape = 1, size = 0.75)
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "EH",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_c("EH of All SNPs using Physical Map")
)

# plot the matrix
png("Results/loci/EH/All_EH_phys.png",
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300)
plots_matrix
dev.off()

###############################################################################
# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_gen_sample_subset")

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(wheat_data$snp$pos, wheat_data$snp$chrom, max)
max_genome_lengths <- c(
  max(chrom_lengths[seq(1, 19, 3)]), # A genome
  max(chrom_lengths[seq(2, 20, 3)]), # B genome
  max(chrom_lengths[seq(3, 21, 3)])  # D genome
)

# create a snp_data frame of all the
wheat_data$snp <- wheat_data$snp %>%
  add_column(eh = calc_eh(wheat_data$genotypes))

haplo_ids <- c(
  haplo_id_1A, haplo_id_2A, haplo_id_4A, haplo_id_5B, haplo_id_6A, haplo_id_6B,
  haplo_id_7A
)
haplo_index_snps <- match(haplo_ids, wheat_data$snp$id)

# create a column in the snp_data set that contains phi values for only
# those markers in the extended haplotypes, NA for all else
wheat_data$snp <- wheat_data$snp %>% mutate(haplo = eh)
wheat_data$snp$haplo[-haplo_index_snps] <- NA

# allows application of same colour to each set of chromosomes
colour_order <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3),
rep(6, 3), rep(7, 3))

# create plots of all markers' eh values on each chromosome
plots <- by(wheat_data$snp, wheat_data$snp$chrom,
  function (data_chrom) {
    chrom_num <- data_chrom$chrom[1]
    data_chrom %>%
      ggplot() +
      xlim(
        0,
        max_genome_lengths[ifelse(chrom_num %% 3, chrom_num %% 3, 3)]
      ) +
      geom_point(
        aes(pos, eh), colour = colours_chroms[colour_order[chrom_num]],
        size = 0.5
      ) +
      geom_point(aes(pos, haplo), shape = 1, size = 0.75)
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "EH",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_c("EH of All SNPs using Genetic Map")
)

# plot the matrix
png("Results/loci/EH/All_EH_gen.png",
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300)
plots_matrix
dev.off()
