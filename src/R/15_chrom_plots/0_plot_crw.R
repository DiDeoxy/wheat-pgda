# load the needed functions and file paths
source("wheat-pgda/src/R/file_paths.R")
library(tidyverse)
library(pgda)
library(GGally)
library(plyr)

wheat_data <- snpgds_parse(phys_gds)

cereba_data <- read_tsv(
  "/workspace/data/intermediate/blast/cereba.txt",
  col_names = c(
    "qseqid", "chrom", "bitscore", "pident", "evalue", "qlen", "length",
    "sstart", "send"
  )
) %>%
  mutate(cereba_pos_mb = (sstart + send) / 2e6) %>%
  select(chrom, cereba_pos_mb) %>%
  rbind.fill(read_csv(file.path(intermediate, "centromeres.csv"))) %>%
  mutate(cent_pos_mb = pos_mb) %>% select(-pos_mb)

# calc the lengths of the different genomes and homoeologous sets
max_phys_lengths <- span_by_chrom(
  wheat_data$snp$chrom, wheat_data$snp$pos
) %>% max_lengths() / 1e6

# create plots of wheat position vs gen pos
plots <- by(cereba_data, cereba_data$chrom,
  function (chrom_data) {
    chrom <- chrom_data$chrom[1]
    plot <- chrom_data %>%
      ggplot() +
      xlim(
        0,
        max_phys_lengths[[
          ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
        ]]
      ) +
      geom_density(aes(cereba_pos_mb), colour = "blue") +
      geom_vline(aes(xintercept = cent_pos_mb), colour = "red")
  }
)

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Position in cM",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_c(
    "Marker Pseudo-Chromosomal vs Genetic Position by\n",
    "Chromosome Coloured by Major Allele Frequency"
  )
)

# plot the matrix
png(
  file.path("results", "wheat_vs_gen_pos_with_cereba.png"),
  family = "Times New Roman", width = 2480, height = 3508
)
plots_matrix
dev.off()