library(tidyverse)
library(plyr)
library(GGally)
library(ggrepel)
library(extrafont)
library(SNPRelate)

# load custom functions
source("src\\R_functions\\funcs_calc_stats.R")
source("src\\R_functions\\colour_sets.R")

# load the data from the gds object
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")

# load the locations of the known genes
known_genes <- read_rds(
  "Data\\Intermediate\\Aligned_genes\\known_genes_groups_corrected.rds")

# function to create a plot for the chromosome of each gene
plot_gene_chrom <- function (gene) {
  amova <- read_rds(str_c("Data\\Intermediate\\Adegenet\\", gene,
                          "_amova.rds"))
  data_phi <- data %>% add_column(phi = phi_markers(amova))
  signif <- quantile(data_phi$phi, prob = 0.975, na.rm = T)
  known_row <- which(known_genes$id == gene)
  data_phi_chrom <- data_phi %>%
                    filter(chrom == known_genes$chrom[known_row])
  ggplot() +
    geom_point(aes(data_phi_chrom$pos, data_phi_chrom$phi),
               colour = colour_set[19], size = 0.5) +
    geom_point(aes(known_genes$pos[known_row], signif),
               shape = 17, colour = colour_set[15]) +
    geom_text_repel(aes(known_genes$pos[known_row], signif,
                        label = known_genes$id[known_row]),
                    angle = 90, hjust = 0, vjust = 1, size = 3,
                    show.legend = FALSE, nudge_y = 0.3) +
    ylim(signif, 1) +
    xlim(0, 650)
}

# making the plots for each gene
plots <- lapply(c("Lr10", "Lr21", "Lr22a", "Lr1", "Lr34"), plot_gene_chrom)

# creating a ggmatrix of the plots
plots_matrix <- ggmatrix(
  plots, nrow = 5, ncol = 1, xlab = "Position in Mb", ylab = "Phi Value",
  yAxisLabels = c("Chr 1A", "Chr 1D", "Chr 2D", "Chr 5D", "Chr 7A"),
  title = str_c("Plots of top 2.5% of Phi Values of Markers for Comparisons ",
                "of Groups Known to Contain Opposite Alleles")
)

# plot the ggmatrix to a png
png("Results\\loci\\amova\\full_amova_genes.png", family = "Times New Roman",
    width = 66, height = 205, pointsize = 5, units = "mm", res = 300)
plots_matrix
dev.off()