library(plyr)
library(tidyverse)
library(GGally)
library(ggrepel)
library(extrafont)

# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/colour_sets.R")
source("src/R_functions/funcs_locus_by_locus.R")

# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

# find the max position of any marker on each genome for xlims
max_genome_lengths <- calc_max_genome_lengths(wheat_data)

# names of groups to be plotted
groups <- c("chrs_csws", "chrs_chrw", "csws_chrw")

# find the phi values of each marker in each Gene and add to data set
wheat_data <- add_group_stat(wheat_data, groups)

# find the extreme threshold for each Gene
extremes <- calc_extremes(wheat_data, groups)

group_extreme_freqs <- calc_group_extreme_freqs(
  wheat_data, extremes
)

# load the gene positions
pheno_genes <- load_groups("pheno_genes.csv") %>%
  mutate(group = "pheno_gene")
resi_genes <- load_groups("resi_genes.csv") %>%
  mutate(group = "resi_gene")

# add the genes in to the data frame
group_extreme_freqs_genes <- group_extreme_freqs %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, group, pos_mb) %>%
  as.tibble()

# read_csv(
#     str_c("Data/Intermediate/Aligned_genes/selected_alignments/pheno_genes.csv"),
#     col_names = c("group", "chrom", "pos", "sleng", "salign", "%id")
#   )

# max(group_extreme_freqs$num, na.rm = T)
# group_extreme_freqs_genes[which(group_extreme_freqs_genes$chrom == "1A"), ]

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
lables <- c(
  "CHRS vs CHRW", "CHRS vs CSWS", "CSWS vs CHRW", "Phenotype Genes",
  "Resistance Genes"
)
legend_title <- "Comparisons and Genes"
plots <- by(
  group_extreme_freqs_genes, group_extreme_freqs_genes$chrom,
  function(chrom_data) {
    chrom <- chrom_data$chrom[1]
    chrom_data %>%
      ggplot() +
      ylim(0, 1) +
      xlim(
        0, 
        max_genome_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ]
      ) +
      geom_point(
        aes(pos_mb, freq, colour = group, size = num), shape = 10
      ) +
      geom_point(
        aes(pos_mb, base, colour = group, shape = group), size = 0.75
      ) +
      geom_text_repel(
        aes(pos_mb, base, colour = group, label = id), angle = 90, hjust = 0,
        vjust = 1, size = 3, segment.colour = "black",
        nudge_y = 0.07,
        nudge_x = ifelse(chrom == "1D", 80,
          ifelse(chrom %in% c("2D", "4A"), -60, 40)
        ),
        show.legend = FALSE
      ) +
      scale_colour_manual(
        legend_title, labels = lables, values = colours_groups_genes,
        limits = levels(as.factor(group_extreme_freqs_genes$group))
      ) +
      scale_shape_manual(
        legend_title, labels = lables, values = c(15, 16, 18, 17, 8),
        limits = levels(as.factor(group_extreme_freqs_genes$group))
      ) +
      labs(colour = "Comparison") +
      scale_size_continuous(
        name = "Number of Markers", trans = "sqrt", limits = c(1, 125)
      )
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = "Average Frequency of Exceptional Jost's D Values in Nearby Markers",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
   title = str_c(
    "Number of Markers in a Region and their Average Frequency of Exceptional ",
    "Jost's D Values"
  ),
  legend = c(1, 1)
)

# plot the matrix
png(str_c("Results/loci/D/comps_D.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()