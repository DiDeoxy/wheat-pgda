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

# find the jost's D values of each marker in each Gene and add to data set
wheat_data <- add_group_stat(wheat_data, groups)

# find the extreme threshold for each Gene
extremes <- calc_extremes(wheat_data, groups)

# create a table of the regions with a high density of extreme markers
group_extreme_freqs <- calc_group_extreme_freqs(
  wheat_data, extremes
)

# load the gene positions
pheno_genes <- load_groups("pheno_genes.csv", base = 0.5) %>%
  mutate(group = "pheno_gene") %>%
  rename(mean_pos_mb = "pos_mb")
resi_genes <- load_groups("resi_genes.csv", base = 0.5) %>%
  mutate(group = "resi_gene") %>%
  rename(mean_pos_mb = "pos_mb")

# add the genes positons to the regions table
group_extreme_freqs_genes <- group_extreme_freqs %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, group, pos_mb) %>%
  as.tibble()

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
lables <- c(
  "CHRS vs CHRW", "CHRS vs CSWS", "CSWS vs CHRW", "Phenotype Genes",
  "Resistance Genes"
)
legend_title <- "Comparisons and Genes"
group_extreme_freqs_genes$mean_pos_mb[is.na(which(group_extreme_freqs_genes$base))]
plots <- by(
  group_extreme_freqs_genes, group_extreme_freqs_genes$chrom,
  function(chrom_data) {
    chrom <- chrom_data$chrom[1]
    chrom_data %>%
      ggplot() +
      ylim(0.5, 1) +
      xlim(
        0, 
        max_genome_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ]
      ) +
      geom_segment(
        aes(
          x = start - 1, y = mean_D, xend = end + 1, yend = mean_D, colour = group,
          size = num_linked
        )
      ) +
      scale_size_continuous(
        "Number of Markers", trans = "sqrt",
        limits = c(1, max(group_extreme_freqs$num_linked, na.rm = T))
      ) +
      geom_point(
        aes(start, mean_D, colour = group), shape = 20, size = 1
      ) +
      geom_point(
        aes(end, mean_D, colour = group), shape = 20, size = 1
      ) +
      geom_point(
        aes(mean_pos_mb, base, colour = group, shape = group), size = 1
      ) +
      geom_text_repel(
        aes(mean_pos_mb, base, colour = group, label = id), angle = 90, hjust = 0,
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
      labs(colour = "Comparison")
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = str_c(
    "Average Jost's D Values in Regions with Nearby Exceptional Markers"
  ),
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Number of Markers in a Region and their Average Jost's D Values",
  legend = c(1, 1)
)

# plot the matrix
png(str_c("Results/loci/D/comps_D.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()

# print our the markers involved in each linked region
comp_gef <- group_extreme_freqs[complete.cases(group_extreme_freqs), ]
genes <- rbind(pheno_genes, resi_genes)
base <- "Results/loci/D/closest_markers"
for (row in 1:nrow(comp_gef)) {
  file_name <- paste(
    cbind(
      comp_gef[row, c("chrom", "mean_pos_mb", "group", "num_linked")] %>%
      round_df(0),
      comp_gef[row, c("freq_extreme", "mean_D")] %>% round_df(2)
    ), collapse = '_'
  )
  linked <- tibble(
    extreme = strsplit(comp_gef[row, "extreme"] %>% as.character(), ' ')[[1]],
    pos_mb = strsplit(comp_gef[row, "pos_mb"] %>% as.character(), ' ')[[1]] %>% as.numeric(),
    Ds = strsplit(comp_gef[row, "Ds"] %>% as.character(), ' ')[[1]],
    ids = strsplit(comp_gef[row, "ids"] %>% as.character(), ' ')[[1]]
  )
  for (row2 in 1:nrow(genes)) {
    if (genes[row2, ]$chrom == comp_gef[row, ]$chrom) {
      print(genes[row2, ]$chrom)
      linked <- linked %>%
        add_row(
          extreme = "NA", pos_mb = genes[row2, ]$mean_pos_mb, Ds = "NA",
          ids = genes[row2, ]$id
        )
    }
  }
  linked <- linked %>% arrange(pos_mb)
  ifelse(! dir.exists(file.path(base)), dir.create(file.path(base)), FALSE)
  write_csv(linked, file.path(base, str_c(file_name, ".csv")))
}

sum(comp_gef$group == "chrs_chrw")
# comp_gef[comp_gef$group == "chrs_chrw", ] %>% print(n = Inf)
sum(comp_gef$group == "chrs_csws")
sum(comp_gef$group == "csws_chrw")
# sum(comp_gef$group == "chrs_chrw" & comp_gef$mean_D > 0.5)
# sum(comp_gef$group == "chrs_chrw" & comp_gef$mean_D > 0.75)
# comp_gef[which(comp_gef$group == "chrs_chrw" & comp_gef$mean_D > 0.75), ] %>% print(n = Inf)