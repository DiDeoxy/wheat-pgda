library(plyr)
library(tidyverse)
library(GGally)
library(ggrepel)
library(extrafont)

# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/funcs_locus_by_locus.R")

# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

# find the max position of any marker on each genome for xlims
max_genome_lengths <- calc_max_genome_lengths(wheat_data)

# names of groups to be plotted
groups <- c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")

# find the phi values of each marker in each Gene and add to data set
wheat_data <- add_group_stat(wheat_data, groups)

# find the extreme threshold for each Gene
extremes <- calc_extremes(wheat_data, groups)

group_extreme_freqs <- calc_group_extreme_freqs(
  wheat_data, extremes
)

# load the gene positions
known_genes <- load_groups("known_genes.csv") %>%
  mutate(group = id)

# add the genes in to the data frame
group_extreme_freqs_known <- group_extreme_freqs %>%
  rbind.fill(known_genes) %>%
  arrange(chrom, group, pos_mb) %>%
  as.tibble()

# max(group_extreme_freqs$num, na.rm = T)

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by Gene or gene type
legend_title <- "Comparisons and Genes"
plots <- by(group_extreme_freqs_known, group_extreme_freqs_known$chrom,
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
        geom_text_repel(
          aes(pos_mb, base, colour = group, label = group), angle = 90, hjust = 0,
          vjust = 1, size = 3, segment.colour = "black", nudge_y = 0.15,
          nudge_x = 60, show.legend = FALSE
        ) +
        labs(colour = "Comparison") +
        scale_size_continuous(
          name = "Number of Markers", trans = "sqrt", limits = c(1, 175)
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
    "Average Frequency and Number of Markers in Regions with Exceptional ",
    "Jost's D values"
  ),
  legend = c(1, 1)
)

plots_matrix[1, 1] <- plots_matrix[1, 1] +
  theme(panel.border = element_rect(fill = NA, colour = "#A3A500", size = 2))
plots_matrix[1, 3] <- plots_matrix[1, 3] +
  theme(panel.border = element_rect(fill = NA, colour = "#00BF7D", size = 2))
plots_matrix[2, 3] <- plots_matrix[2, 3] +
  theme(panel.border = element_rect(fill = NA, colour = "#00B0F6", size = 2))
plots_matrix[5, 3] <- plots_matrix[5, 3] +
  theme(panel.border = element_rect(fill = NA, colour = "#F8766D", size = 2))
plots_matrix[7, 3] <- plots_matrix[7, 3] +
  theme(panel.border = element_rect(fill = NA, colour = "#E76BF3", size = 2))

# plot the matrix
png(str_c("Results/loci/D/known_genes_D_freq_point.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()