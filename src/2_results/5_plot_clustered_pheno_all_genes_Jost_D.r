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
groups <- c("chrs_chrw", "chrs_csws", "csws_chrw")
# groups <- c("chrs_chrw_csws")

# find the jost's D values of each marker in each Gene and add to data set
wheat_data <- add_group_stat(wheat_data, groups)
wheat_data$snp <- wheat_data$snp %>%
  gather(group, D, groups)

# load the gene positions
pheno_genes <- load_groups("pheno_genes.csv", base = 1) %>%
  mutate(group = "pheno_gene")
resi_genes <- load_groups("resi_genes.csv", base = 1) %>%
  mutate(group = "resi_gene")

# # add the genes positons to the regions table
wheat_data$snp <- wheat_data$snp %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, group, pos_mb) %>%
  as_tibble()

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
lables <- c(
  "CHRS vs CHRW", "CHRS vs CSWS", "CSWS vs CHRW", "Phenotype Genes",
  "Resistance Genes"
)
legend_title <- "Comparisons and Genes"

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
        aes(pos_mb, D, colour = group), shape = 16, size = 1, alpha = 1/8
      ) +
      geom_smooth(
        aes(pos_mb, D, colour = group), method = "loess", span = 0.06, size = 0.5
      ) +
      geom_point(
        aes(pos_mb, base, colour = group, shape = group), size = 1
      ) +
      geom_text_repel(
        aes(pos_mb, base, colour = group, label = id), angle = 90, hjust = 0,
        vjust = -1, size = 2, fontface = "bold",
        nudge_y = -0.07,
        nudge_x = ifelse(chrom == "1D", 80,
          ifelse(chrom %in% c("2D", "4A"), -60, 40)
        ),
        show.legend = FALSE
      ) +
      scale_colour_manual(
        legend_title, labels = lables, values = colours_groups_genes,
        limits = levels(as.factor(wheat_data$snp$group))
      ) +
      scale_shape_manual(
        legend_title, labels = lables, values = c(15, 16, 18, 17, 8),
        limits = levels(as.factor(wheat_data$snp$group))
      )
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
png("Results/loci/D/comps_D_test.png",
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()

# ################################################################################
# # printing out the data in various ways
genes <- rbind(pheno_genes, resi_genes) %>% as_tibble()
base <- "Results/loci/D/closest_markers"

wheat_data$snp <- wheat_data$snp %>%
  add_column(chrs_chrw_exceptional_Ds = (wheat_data$snp$chrs_chrw >= extremes["chrs_chrw"])) %>%
  add_column(chrs_csws_exceptional_Ds = (wheat_data$snp$chrs_csws >= extremes["chrs_csws"])) %>%
  add_column(csws_chrw_exceptional_Ds = (wheat_data$snp$csws_chrw >= extremes["csws_chrw"]))
wheat_data$snp <- wheat_data$snp[, c(1:4, 8:10, 5:7)]
# print out all markers on each chromosome
ifelse(! dir.exists(file.path(base)), dir.create(file.path(base)), FALSE)
blah <- by(wheat_data$snp, wheat_data$snp$chrom, function (chrom) {
  for (i in 1:nrow(genes)) {
    if (genes$chrom[i] == chrom$chrom[1]) {
      chrom <- chrom %>%
        add_row(
          id = genes$id[i], chrom = genes$chrom[i], pos_mb = genes$mean_pos_mb[i]
        )
    }
  }
  write_csv(chrom %>% arrange(pos_mb), file.path(base, str_c(chrom$chrom[1], "_all_markers_josts_D.csv")))
})

###########
comp_gef <- group_extreme_freqs[complete.cases(group_extreme_freqs), ]

# print our summaries of the regions identified in each comp with genes
ifelse(! dir.exists(file.path(base)), dir.create(file.path(base)), FALSE)
blah <- by(comp_gef, comp_gef$chrom, function (chrom) {
  for (row in 1:nrow(genes)) {
    if (genes[row, ]$chrom == chrom$chrom[1]) {
      chrom <- chrom %>%
        add_row(
          chrom = genes[row, ]$chrom, start = genes[row, ]$mean_pos_mb, group = genes[row, ]$id
        )
    }
  }
  write_csv(chrom[, 1:7] %>% arrange(start, group), file.path(base, str_c(chrom$chrom[1], "_regions_summary.csv")))
})

# # print out for each chromosome the markers of each region with line gaps for
# # when overlaps occur or change or a new region starts
# print_ovlps_by_chrom(all_ovlps_markers(comp_gef))

# # identify which other groups overlap each region of a specific group
# comp_ovlps(comp_gef, "chrs_chrw") %>% unlist() %>% table()
# comp_ovlps(comp_gef, "chrs_csws") %>% unlist() %>% table()
# comp_ovlps(comp_gef, "csws_chrw") %>% unlist() %>% table()

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
