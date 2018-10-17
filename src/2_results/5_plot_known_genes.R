library(tidyverse)
library(plyr)
library(GGally)
library(ggrepel)
library(extrafont)
library(SNPRelate)

# load custom functions
source("src\\R_functions\\funcs_plot_loci.R")
source("src\\R_functions\\funcs_calc_map_stats.R")

# load the data from the gds object
wheat_data <- parse_gds("phys_subset_sample")

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(wheat_data$snp$pos_mb, wheat_data$snp$chrom, max)
max_genome_lengths <- data.frame(A = max(chrom_lengths[seq(1, 19, 3)]),
                                 B = max(chrom_lengths[seq(2, 20, 3)]),
                                 D = max(chrom_lengths[seq(3, 21, 3)]))

# # find the phi values of each marker in each Gene and add to data set
# for (group in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
#   amova <- read_rds(str_c("Data\\Intermediate\\Adegenet\\", group,
#                           "_amova.rds"))
#   wheat_data$snp <- wheat_data$snp %>%
#     add_column(!!group := retrieve_phi(amova))
# }

# find the phi values of each marker in each Gene and add to data set
for (group in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
  gene_diff_stats <- read_rds(str_c("Data\\Intermediate\\mmod\\", group,
    "_diff_stats.rds"))
  gene_diff_stats <- gene_diff_stats[[1]] %>% as.tibble()
  wheat_data$snp <- wheat_data$snp %>%
    add_column(!!group := gene_diff_stats$D)
}

# find the significance threshold for each Gene
signifs <- wheat_data$snp %>%
  summarise(
    Lr34 = quantile(Lr34, prob = 0.975, na.rm = TRUE),
    Lr22a = quantile(Lr22a, prob = 0.975, na.rm = TRUE),
    Lr21 = quantile(Lr21, prob = 0.975, na.rm = TRUE),
    Lr10 = quantile(Lr10, prob = 0.975, na.rm = TRUE),
    Lr1 = quantile(Lr1, prob = 0.975, na.rm = TRUE)
  )

# turn each value less than the significance threshold for each Gene to
# NA
wheat_data$snp$Lr34[wheat_data$snp$Lr34 < signifs$Lr34] <- NA
wheat_data$snp$Lr22a[wheat_data$snp$Lr22a < signifs$Lr22a] <- NA
wheat_data$snp$Lr21[wheat_data$snp$Lr21 < signifs$Lr21] <- NA
wheat_data$snp$Lr10[wheat_data$snp$Lr10 < signifs$Lr10] <- NA
wheat_data$snp$Lr1[wheat_data$snp$Lr1 < signifs$Lr1] <- NA

# make the data set long for easier plotting
wheat_data$snp_long <- wheat_data$snp %>%
  gather(Gene, phi, c(Lr34, Lr22a, Lr21, Lr10, Lr1))

# load the gene positions
known_genes <- read_csv(
  "Data\\Intermediate\\Aligned_genes\\selected_alignments\\known_genes.csv",
  col_names = c("Gene", "chrom", "pos", "sleng", "salign", "%id")
) %>%
  select(Gene, chrom, pos) %>%
  mutate(pos_mb = pos / 1000000)

# convert known genes chromosome  names to appropriate integer
chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()
for (i in seq_along(chroms)) {
  known_genes$chrom[known_genes$chrom == chroms[i]] <- i
}
known_genes$chrom <- as.integer(known_genes$chrom)

# add min phi values of plotting to each gene so they appear at bottom of plots
known_genes <- cbind(known_genes, min_signif = min(signifs))
known_genes <- cbind(known_genes, gene_signif = 0)
for (i in 1:nrow(known_genes)) {
  known_genes$gene_signif[i] <- signifs[, known_genes$Gene[i]]
}

# add the genes in to the data frame
wheat_data$snp_long_genes <- wheat_data$snp_long %>%
  rbind.fill(known_genes) %>%
  arrange(chrom, Gene, pos) %>%
  as.tibble()

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by Gene or gene type
lables <- c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")
legend_title <- "Comparisons and Genes"
plots <- by(
  wheat_data$snp_long_genes, wheat_data$snp_long_genes$chrom,
  function(data_chrom) {
    chrom_num <- as.integer(data_chrom$chrom[1])
    data_chrom %>%
      ggplot() +
      ylim(min(signifs), 1) +
      xlim(
        0,
        max_genome_lengths[, ifelse(chrom_num %% 3, chrom_num %% 3, 3)]
      ) +
      geom_point(
        aes(pos_mb, phi, colour = Gene, shape = Gene), size = 0.5
      ) +
      geom_text_repel(aes(pos_mb, min_signif, colour = Gene,
        label = Gene), angle = 90, hjust = 0, vjust = 1, size = 3,
        segment.colour = "black", nudge_y = 0.3, nudge_x = 40,
        show.legend = FALSE
      )
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Phi",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Markers in Top 2.5% of AMOVA Phi Values For Resistance Gene Groups",
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
png(str_c("Results\\loci\\amova\\known_genes.png"),
  family = "Times New Roman", width = 210, height = 277, pointsize = 5,
  units = "mm", res = 300
)
print(plots_matrix + theme(legend.position = "bottom"))
dev.off()

###############################################################################
# find markers closest to homoeologs responsible for resistance phenotype
wheat_data <- parse_gds("phys_subset_sample")

# for (group in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
#   amova <- read_rds(str_c("Data\\Intermediate\\Adegenet\\", group,
#                           "_amova.rds"))
#   wheat_data$snp <- wheat_data$snp %>%
#     add_column(!!group := retrieve_phi(amova))
# }

# find the phi values of each marker in each Gene and add to data set
# for (group in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
#   gene_genind <- read_rds(str_c("Data\\Intermediate\\Adegenet\\", group,
#     "_genind.rds"))
#   stats1 <- diff_stats(gene_genind)[[1]] %>% as.tibble()
#   wheat_data$snp <- wheat_data$snp %>%
#     add_column(!!group := stats1$Gprime_st)
# }

known_genes_responsible <- cbind(
  known_genes[c(1, 6, 9, 11, 15), ], dist = c(20, 1.92, 7.4, 0.193, 11.25)
) %>% as.tibble()

for (i in 1:nrow(known_genes_responsible)) {
  row <- known_genes_responsible[i, ]
  nearby_markers <- wheat_data$snp %>%
    filter(chrom == row$chrom & abs(pos_mb - row$pos_mb) < row$dist) %>%
    transmute(Marker = id, Chrom = chrom, Distance = pos_mb - row$pos_mb,
      phi := !!as.symbol(row$Gene))
  nearby_signif_markers <- wheat_data$snp %>%
    filter(chrom == row$chrom & abs(pos_mb - row$pos_mb) < row$dist &
      !!as.symbol(row$Gene) > row$gene_signif) %>%
    transmute(Marker = id, Chrom = chrom, Distance = pos_mb - row$pos_mb,
      phi := !!as.symbol(row$Gene))
  write_csv(nearby_markers,
    str_c("Results\\closest_markers\\", row$Gene, "_full.csv"))
  write_csv(nearby_signif_markers,
    str_c("Results\\closest_markers\\", row$Gene, "_signif.csv"))
}

# Lr 10
(0.310631 + 0.327516 + 0.310631 + 0.310631 + 0.293038 + 0.293038 + 0.347754 +
0.351068 + 0.275815 + 0.310631 + 0.346957 + 0.346957 + 0.310631 + 0.398613 + 
0.39781) / 15

# Lr21
(0.340514 + 0.374857 + 0.374857 + 0.940899 + 0.62186 + 0.62186 + 0.62186 +
0.594675 + 0.62186) / 9