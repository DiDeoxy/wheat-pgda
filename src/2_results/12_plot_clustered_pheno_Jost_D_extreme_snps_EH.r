library(tidyverse)
library(plyr)
library(GGally)
# install.packages("gridExtra")
library(gridExtra)
library(SNPRelate)
library(extrafont)

# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/funcs_locus_by_locus.R")
source("src/R_functions/colour_sets.R")

# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

# find the max position of any marker on each genome for xlims
max_genome_lengths <- calc_max_genome_lengths(wheat_data)

# names of groups to be plotted
groups <- c("chrs_chrw", "chrs_csws", "csws_chrw")

# load the hdbscan clusters
cluster <- read_rds("Data/Intermediate/hdbscan/wheat_hdbscan.rds")$cluster

# find the indexes of the three phenotypes in the three major clusters
sample_indices <- list(
  chrs = which(wheat_data$sample$annot$pheno == "HRS" & cluster == 5),
  chrw = which(wheat_data$sample$annot$pheno == "HRW" & cluster == 1),
  csws = which(wheat_data$sample$annot$pheno == "SWS" & cluster == 2)
)

# find the extreme threshold for each gene
# find the phi values of each marker in each Gene and add to data set
wheat_data <- add_group_stat(wheat_data, groups)
extremes <- calc_extremes(wheat_data, groups)

# find the markers which are extreme in each comparison
marker_indices <- list(
  chrs_chrw = which(wheat_data$snp$chrs_chrw >= extremes[["chrs_chrw"]]),
  chrs_csws = which(wheat_data$snp$chrs_csws >= extremes[["chrs_csws"]]),
  csws_chrw = which(wheat_data$snp$csws_chrw >= extremes[["csws_chrw"]])
)

# calc eh for the extreme markers in each comparison for each group in that
# comparison
eh <- list(
  chrs_chrw = tibble(
    chrs = calc_eh(
      wheat_data$genotypes[marker_indices$chrs_chrw, sample_indices$chrs]
    ),
    chrw = calc_eh(
      wheat_data$genotypes[marker_indices$chrs_chrw, sample_indices$chrw]
    )
  ),
  chrs_csws = tibble(
    chrs = calc_eh(
      wheat_data$genotypes[marker_indices$chrs_csws, sample_indices$chrs]
    ),
    csws = calc_eh(
      wheat_data$genotypes[marker_indices$chrs_csws, sample_indices$csws]
    )
  ),
  csws_chrw = tibble(
    csws = calc_eh(
      wheat_data$genotypes[marker_indices$csws_chrw, sample_indices$csws]
    ),
    chrw = calc_eh(
      wheat_data$genotypes[marker_indices$csws_chrw, sample_indices$chrw]
    )
  )
)

lim <- c(0, 0.3)
# plot the EH values for each group in each comparison on each chromosome
plots <- list(
  chrs_chrw = ggplot() +
    xlim(lim) +
    ylim(lim) +
    geom_density_2d(
      aes(eh$chrs_chrw$chrs, eh$chrs_chrw$chrw),
      colour = colours_groups_genes[1]) + 
    labs(
      title = "EHs of Markers With\nExceptional Jost's D\nin CHRS vs CHRW",
      x = "EH CHRS", y = "EH CHRW"
    ),
  chrs_csws = ggplot() +
    xlim(lim) +
    ylim(lim) +
    geom_density_2d(
      aes(eh$chrs_csws$chrs, eh$chrs_csws$csws),
      colour = colours_groups_genes[2]) +
    labs(
      title = "EHs of Markers With\nExceptional Jost's D\nin CHRS vs CSWS",
      x = "EH CHRS", y = "EH CSWS"
    ),
  csws_chrw = ggplot() +
    xlim(lim) +
    ylim(lim) +
    geom_density_2d(
      aes(eh$csws_chrw$csws, eh$csws_chrw$chrw),
      colour = colours_groups_genes[3]) +
    labs(
      title = "EHs of Markers With\nExceptional Jost's D\nin CSWS vs CHRW",
      x = "EH CSWS", y = "EH CHRW"
    )
)

# plot
png(str_c("Results/loci/EH/extreme_Jost_D_eh_comps.png"),
    family = "Times New Roman", width = 210, height = 75, pointsize = 3,
    units = "mm", res = 300)
grid.arrange(grobs = plots, nrow = 1)
dev.off()