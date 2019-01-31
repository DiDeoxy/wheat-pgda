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
groups <- c("chrs_chrw_csws")

# find the jost's D values of each marker in each Gene and add to data set
wheat_data$snp <- add_group_stat(wheat_data$snp, groups) %>% as.tibble()
wheat_data$snp <- wheat_data$snp %>%
  gather(group, D, groups)

# load the gene positions
pheno_genes <- load_groups("pheno_genes.csv", base = 1) %>%
  mutate(type = "pheno_gene")
resi_genes <- load_groups("resi_genes.csv", base = 1) %>%
  mutate(type = "resi_gene")

# # add the genes positons to the regions table
wheat_data$snp <- wheat_data$snp %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, group, pos_mb) %>%
  as_tibble()

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
lables <- c("Phenotype Genes", "Resistance Genes")

legend_title <- "Genes"
chroms_order <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>% 
    t() %>% as.vector()
colour_order <- rbind(1:7, 1:7, 1:7) %>% as.vector()

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
        aes(pos_mb, D),
        colour = colours_chroms[colour_order[which(chroms_order == chrom)]],
        shape = 16, size = 1, alpha = 1/4
      ) +
      geom_smooth(
        aes(pos_mb, D),
        colour = colour_set[22],
        method = "loess", span = 0.06, size = 0.375,
      ) +
      geom_point(
        aes(pos_mb, base, colour = type, shape = type), size = 1
      ) +
      geom_text_repel(
        aes(pos_mb, base, colour = type, label = id), angle = 90, hjust = 0,
        vjust = -1, size = 3, fontface = "bold",
        nudge_y = -0.07,
        nudge_x = ifelse(chrom == "1D", 80,
          ifelse(chrom %in% c("2D", "4A"), -60, 40)
        ),
        show.legend = FALSE
      ) +
      scale_colour_manual(
        legend_title, labels = lables, values = colours_groups_genes[4:5],
        limits = levels(as.factor(wheat_data$snp$type))
      ) +
      scale_shape_manual(
        legend_title, labels = lables, values = c(15, 16, 18, 17, 8)[4:5],
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

# ################################################################################
# # printing out the data in various ways
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")
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
# }) %>% t() %>% as.tibble()
# colnames(major_allele_freq_pops) <- c("chrs", "chrw", "csws")

# load the gene positions
pheno_genes <- load_groups("pheno_genes.csv", base = 1) %>%
  select(-"base")
resi_genes <- load_groups("resi_genes.csv", base = 1) %>%
  select(-"base")

wheat_data$snp <- wheat_data$snp %>%
  add_group_stat("chrs_chrw_csws") %>%
  gather(group, D, "chrs_chrw_csws") %>%
  select(-c(group, pos)) %>%
  mutate(D = D %>% round(4)) %>%
  cbind(major_allele_freq_pops %>% round(4)) %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, pos_mb) %>%
  as_tibble()

# print out all markers on each chromosome
base <- "Results/loci/D/by_chrom"
ifelse(! dir.exists(file.path(base)), dir.create(file.path(base)), FALSE)
blah <- by(wheat_data$snp, wheat_data$snp$chrom, function (chrom) {
  write_csv(chrom, file.path(base, str_c(chrom$chrom[1], "_all_markers_Ds_and_major_allele_freqs_by_pop_with_genes.csv")))
})

# ################################################################################
# comp_gef %>% print(n = Inf)

# sum(comp_gef$group == "chrs_chrw")
# sum(comp_gef$group == "chrs_csws")
# sum(comp_gef$group == "csws_chrw")
# comp_gef[which(comp_gef$group == "chrs_chrw"), ] %>% print(n = Inf)

# # markers near Lr22A
# D2 <- which(wheat_data$snp$chrom == "2D")
# sum(wheat_data$snp$pos_mb[D2] >= 1.669121 & wheat_data$snp$pos_mb[D2] <= 31.04408)
# sum(wheat_data$snp$pos_mb[D2] >= 31.887715 & wheat_data$snp$pos_mb[D2] <= 36)


# A3 <- which(wheat_data$snp$chrom == "3A")
# wheat_data$snp[A3, ][
#   wheat_data$snp$pos_mb[A3] >= (100.881838 - 15)
#   & wheat_data$snp$pos_mb[A3] <= (100.881838 + 15),
# ]

# # markers near Lr34 4A
# B3 <- which(wheat_data$snp$chrom == "3B")
# sum(
#   wheat_data$snp$pos_mb[B3] >= 0
#   & wheat_data$snp$pos_mb[B3] <= 17
# )
# wheat_data$snp[B3, ][
#   wheat_data$snp$pos_mb[B3] >= (757.919031 - 6)
#   & wheat_data$snp$pos_mb[B3] <= (757.919031 + 6),
# ]


# # markers near Tamyb10-D1
# D3 <- which(wheat_data$snp$chrom == "3D")
# wheat_data$snp[D3, ][
#   wheat_data$snp$pos_mb[D3] >= (570.80153 - 10)
#   & wheat_data$snp$pos_mb[D3] <= (570.80153 + 10), 
# ]


# # markers near Lr34 4A
# A4 <- which(wheat_data$snp$chrom == "4A")
# sum(wheat_data$snp$pos_mb[A4] >= 625.862234 & wheat_data$snp$pos_mb[A4] <= 744.309803)]

# # markers near Lr34 4A
# D4 <- which(wheat_data$snp$chrom == "4D")
# sum(
#   wheat_data$snp$pos_mb[D4] >= 18.781996
#   & wheat_data$snp$pos_mb[D4] <= ((20.579698 - 18.781996) + 18.781996)
# )
# wheat_data$snp[D4, ][
#   wheat_data$snp$pos_mb[D4] >= (20.579698 - 6)
#   & wheat_data$snp$pos_mb[D4] <= (20.579698 + 6),
# ]

# B6 <- which(wheat_data$snp$chrom == "6B")
# sum(
#   wheat_data$snp$pos_mb[B6] >= 701
#   & wheat_data$snp$pos_mb[B6] <= 731
# )

# # markers near Lr34 4A
# B7 <- which(wheat_data$snp$chrom == "7B")
# sum(
#   wheat_data$snp$pos_mb[B7] >= 0
#   & wheat_data$snp$pos_mb[B7] <= 17
# )
# wheat_data$snp[B7, ][
#   wheat_data$snp$pos_mb[B7] >= (9.70259 - (15.592659 - 9.70259))
#   & wheat_data$snp$pos_mb[B7] <= (15.592659),
# ]

# # chrs vs CHRW
# library(adegenet)
# comp_genind <- read_rds("Data/Intermediate/Adegenet/chrs_chrw_genind.rds")
# markers <- c("wsnp_Ex_c26128_35374652", "BS00074287_51", "BS00042456_51", "Kukri_c3338_271")
# for (marker in markers) {
#   print(marker)
#   tibble(
#     pop = comp_genind@pop, allele = comp_genind@tab[, str_c(marker, ".A")]
#   ) %>% table() %>% print()
# }

# library(adegenet)
# comp_genind <- read_rds("Data/Intermediate/Adegenet/chrs_csws_genind.rds")
# markers <- c("wsnp_Ex_rep_c68169_66940235")
# for (marker in markers) {
#   print(marker)
#   tibble(
#     pop = comp_genind@pop, allele = comp_genind@tab[, str_c(marker, ".A")]
#   ) %>% table() %>% print()
# }
