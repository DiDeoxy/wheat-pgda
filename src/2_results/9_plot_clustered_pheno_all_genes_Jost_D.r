library(plyr)
library(tidyverse)
library(GGally)
library(ggrepel)
library(extrafont)

# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/colour_sets.R")

# load the data from the gds object
wheat_data <- parse_gds("phys_subset_sample")

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(wheat_data$snp$pos_mb, wheat_data$snp$chrom, max)
max_genome_lengths <- c(
  A = max(chrom_lengths[seq(1, 19, 3)]),
  B = max(chrom_lengths[seq(2, 20, 3)]),
  D = max(chrom_lengths[seq(3, 21, 3)])
)

# add columns of the D values of each comparison
for (group in c("chrs_csws", "chrs_chrw", "csws_chrw")) {
  comp_Jost_D <- read_rds(
    str_c("Data/Intermediate/mmod/", group, "_Jost_D.rds")
  )[[1]]
  wheat_data$snp <- wheat_data$snp %>% add_column(!!group := comp_Jost_D)
}

# find the top 2.5% quantile for each of the three comparisons
extremes <- wheat_data$snp %>%
  summarise(
    chrs_csws = quantile(chrs_csws, prob = 0.975, na.rm = T),
    chrs_chrw = quantile(chrs_chrw, prob = 0.975, na.rm = T),
    csws_chrw = quantile(csws_chrw, prob = 0.975, na.rm = T)
  )

# turn each value less than the top 2.5% quantile for the comparison to NA
wheat_data$snp$chrs_csws[
  which(wheat_data$snp$chrs_csws < extremes$chrs_csws)] <- NA
wheat_data$snp$chrs_chrw[
  which(wheat_data$snp$chrs_chrw < extremes$chrs_chrw)] <- NA
wheat_data$snp$csws_chrw[
  which(wheat_data$snp$csws_chrw < extremes$csws_chrw)] <- NA

# make the data set long for easier plotting
wheat_data$snp_long <- wheat_data$snp %>%
  gather(comparison, D, c(chrs_csws, chrs_chrw, csws_chrw))

# load gene data
pheno_genes <- read_csv(
  "Data/Intermediate/Aligned_genes/selected_alignments/pheno_genes.csv",
  col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
) %>%
  select(id, chrom, pos) %>%
  mutate(pos_mb = pos / 1e6, comparison = "pheno_gene")
resi_genes <- read_csv(
  "Data/Intermediate/Aligned_genes/selected_alignments/resi_genes.csv",
  col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
) %>%
  select(id, chrom, pos) %>%
  mutate(pos_mb = pos / 1e6, comparison = "resi_gene")
# join them together
genes <- pheno_genes %>% full_join(resi_genes)

# convert known genes chromosome  names to appropriate integer
chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()
for (i in seq_along(chroms)) {
  genes$chrom[genes$chrom == chroms[i]] <- i
}
genes$chrom <- as.integer(genes$chrom)

# add min phi values of plotting to each gene so they appear at bottom of plots
genes <- cbind(genes, min_extreme = min(extremes))

# add the genes in and make longer
wheat_data$snp_long_genes <- wheat_data$snp_long %>%
  full_join(genes) %>%
  arrange(chrom, comparison, pos_mb)

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
lables <- c(
  "CHRS vs CHRW", "CHRS vs CSWS", "CSWS vs CHRW", "Phenotype Genes",
  "Resistance Genes"
)
legend_title <- "Comparisons and Genes"
plots <- by(
  wheat_data$snp_long_genes, wheat_data$snp_long_genes$chrom,
  function(chrom) {
    chrom_num <- chrom$chrom[1]
    chrom %>%
      ggplot() +
      ylim(min(extremes), 1) +
      xlim(
        0,
        max_genome_lengths[ifelse(chrom_num %% 3, chrom_num %% 3, 3)]
      ) +
      geom_point(
        aes(pos_mb, D, colour = comparison, shape = comparison), size = 0.75
      ) +
      geom_point(aes(pos_mb, min_extreme, colour = comparison,
        shape = comparison), size = 0.75
      ) +
      geom_text_repel(
        aes(pos_mb, min_extreme, colour = comparison, label = id),
        angle = 90, hjust = 0, vjust = 1, size = 3, segment.colour = "black",
        nudge_y = 0.07, show.legend = FALSE,
        nudge_x =  ifelse(chrom_num == 3, 80,
          ifelse(chrom_num %in% c(6, 10), -60, 40)
        )
      ) +
      scale_colour_manual(
        legend_title, labels = lables, values = colours_comparisons_genes,
        limits = levels(as.factor(wheat_data$snp_long_genes$comparison))
      ) +
      scale_shape_manual(
        legend_title, labels = lables, values = c(15, 16, 18, 17, 8),
        limits = levels(as.factor(wheat_data$snp_long_genes$comparison))
      )
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Jost D",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Markers in Top 2.5% of Jost D Values By Comparison",
  legend = c(1, 1)
)

# plot the matrix
png(str_c("Results/loci/D/comps_D.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom")
dev.off()