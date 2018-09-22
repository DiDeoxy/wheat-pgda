library(tidyverse)
# library(plyr)
# library(GGally)
# library(ggrepel)
# library(extrafont)
library(SNPRelate)

# # load custom functions
# source("src\\R_functions\\funcs_calc_stats.R")
# source("src\\R_functions\\colour_sets.R")

# # load the data from the gds object
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")
table(mc)

# # load the locations of the known genes
# known_genes <- read_rds(
#   "Data\\Intermediate\\Aligned_genes\\known_genes_groups_corrected.rds")

# find the markers closest to the gene in full and > signif set
closest_markers <- function(gene, dist) {
  amova <- read_rds(str_c(
    "Data\\Intermediate\\Adegenet\\", gene,
    "_amova.rds"
  ))
  data_phi <- data %>% add_column(phi = phi_markers(amova))
  signif <- quantile(data_phi$phi, prob = 0.975, na.rm = T)
  gene_row <- which(known_genes$id == gene)
  gene_chrom <- known_genes$chrom[gene_row]
  gene_pos <- known_genes$pos[gene_row]
  full <- data_phi %>%
    filter(chrom == gene_chrom &
      abs(pos - gene_pos) < dist) %>%
    transmute(id, chrom, pos - gene_pos, phi)
  signif <- data_phi %>%
    filter(chrom == gene_chrom & abs(pos - gene_pos) < dist &
      phi > signif) %>%
    transmute(id, chrom, pos - gene_pos, phi)
  write_csv(full, str_c("Results\\closest_markers\\", gene, "_full.csv"))
  write_csv(signif, str_c("Results\\closest_markers\\", gene, "_signif.csv"))
}

# finds the markers
mapply(
  closest_markers, c("Lr10", "Lr21", "Lr22a", "Lr1", "Lr34"),
  c(0.5, 2, 5, 0.2, 12)
)

# # function to create a plot for the chromosome of each gene
# plot_gene_chrom <- function (gene) {
#   amova <- read_rds(str_c("Data\\Intermediate\\Adegenet\\", gene,
#     "_amova.rds"))
#   data_phi <- data %>% add_column(phi = phi_markers(amova))
#   signif <- quantile(data_phi$phi, prob = 0.975, na.rm = T)
#   gene_row <- which(known_genes$id == gene)
#   gene_chrom <- known_genes$chrom[gene_row]
#   data_phi_chrom_signif <- data_phi %>%
#     filter(chrom == gene_chrom & phi > signif)
#   ggplot() +
#     geom_point(aes(data_phi_chrom_signif$pos, data_phi_chrom_signif$phi),
#       colour = colour_set[19], size = 0.5) +
#     geom_point(aes(known_genes$pos[gene_row], 0.2),
#       shape = 17, colour = colour_set[15]) +
#     ylim(0.2, 1) +
#     xlim(0, 650)
# }

# # making the plots for each gene
# plots <- lapply(c("Lr10", "Lr21", "Lr22a", "Lr1", "Lr34"), plot_gene_chrom)

# # creating a ggmatrix of the plots
# plots_matrix <- ggmatrix(
#   plots, nrow = 1, ncol = 5, xlab = "Position in Mb", ylab = "Phi Value",
#   xAxisLabels = c("Chr 1A: Lr10", "Chr 1D: Lr21", "Chr 2D: Lr22a",
#     "Chr 5D: Lr 1", "Chr 7A: Lr34"),
#   title = str_c("Figure 5: Top 2.5% of Phi Statistics by Comparison")
# )

# # plot the ggmatrix to a png
# png("Results\\loci\\amova\\full_amova_genes.png", family = "Times New Roman",
#     width = 200, height = 55, pointsize = 5, units = "mm", res = 300)
# plots_matrix
# dev.off()