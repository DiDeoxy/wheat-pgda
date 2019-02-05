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

cluster <- read_rds("Data/Intermediate/hdbscan/wheat_hdbscan.rds")$cluster

index_chrs_chrw_csws <- which(
  (wheat_data$sample$annot$pheno == "HRS" & cluster == 5)
  | (wheat_data$sample$annot$pheno == "HRW" & cluster == 1)
  | (wheat_data$sample$annot$pheno == "SWS" & cluster == 2)
)
cpheno <- wheat_data$sample$annot$pheno[index_chrs_chrw_csws] %>% factor()
geno <- wheat_data$genotype[, index_chrs_chrw_csws]

# major_allele_freq_pops <- apply(geno, 1, function (marker) {
#   total_geno <- table(marker)
#   major <- which.max(total_geno) %>% names()
#   by(marker, cpheno, function (pop) {
#     pop_geno <- table(pop)
#     if (all(c("0", "2") %in% names(pop_geno))) {
#       pop_geno[[major]] / sum(pop_geno[["0"]], pop_geno[["2"]])
#     } else if (major %in% names(pop_geno)) {
#       1
#     } else {
#       0
#     }
#   }) %>% rbind()
# }) %>% t() %>% as_tibble()
# colnames(major_allele_freq_pops) <- c("chrs", "chrw", "csws")

# load the gene positions
pheno_genes <- load_groups("pheno_genes.csv", base = 1) %>%
  mutate(gene_type = "Phenotype Genes")
resi_genes <- load_groups("resi_genes.csv", base = 1) %>%
  mutate(gene_type = "Resistance Genes")

group <- "chrs_chrw_csws"
# add the genes positons to the regions table
wheat_data$snp <- wheat_data$snp %>%
  add_group_stat(group) %>%
  gather(group, D, group) %>%
  cbind(major_allele_freq_pops %>% round(4)) %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, group, pos_mb) %>%
  rowwise() %>%
  mutate(comp_type = 
    ifelse(chrw >= 0.5 && csws >= 0.5, "None",
      ifelse(chrw >= 0.5 && csws <= 0.5, "CSWS vs CHRS & CHRW",
        ifelse(chrw <= 0.5 && csws >= 0.5, "CHRW vs CHRS & CSWS", 
          ifelse(chrw <= 0.5 && csws <= 0.5, "CHRS vs CHRW & CSWS", NA)
        )
      )
    )
  ) %>%
  mutate(type = pmin(gene_type, comp_type, na.rm = TRUE)) %>%
  select(-c(gene_type, comp_type))

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
legend_title <- "Marker Types & Genes"
plots <- by(
  wheat_data$snp, wheat_data$snp$chrom, function(chrom_groups) {
    chrom <- chrom_groups$chrom[1]
    chrom_groups %>%
      ggplot() +
      ylim(0, 1) +
      xlim(
        0, 
        max_genome_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ]
      ) +
      geom_point(
        aes(pos_mb, D, colour = type),
        shape = 16, size = 1.25, alpha = 1/2
      ) +
      geom_smooth(
        aes(pos_mb, D),
        colour = colour_set[22],
        method = "loess", span = 0.06, size = 0.375,
      ) +
      geom_point(
        aes(pos_mb, base, colour = type), size = 1
      ) +
      geom_text_repel(
        aes(pos_mb, base, colour = type, label = id), angle = 90, hjust = 0,
        vjust = -1, size = 3, fontface = "bold",
        nudge_y = -0.07,
        nudge_x =
          ifelse(chrom %in% c("1A", "1D"), 80,
            ifelse(chrom %in% c("2D", "4A"), -60, 40)
          ),
        show.legend = FALSE
      ) + 
      scale_colour_manual(
        legend_title, values = colours_comps_genes,
        limits = levels(as.factor(wheat_data$snp$type))
      )
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = str_c(
    "Normalized Jost's D Value"
  ),
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Normalized Jost's D Values with Loess Curve By Chromosome ",
  legend = c(1, 1)
)

# plot the matrix
png("Results/loci/D/comps_D_test.png",
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()