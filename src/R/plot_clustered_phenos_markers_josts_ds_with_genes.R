# load file paths, colours, and functions
source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "arrange", "mutate", "rowwise", "select")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_point", "geom_smooth", "guide_legend",
  "scale_colour_manual", "theme", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(pgda, "load_genes", "snpgds_parse")
import::from(plyr, "rbind.fill")
import::from(readr, "read_rds", "write_csv")
import::from(stringr, "str_c")
import::from(tibble, "add_column", "add_row", "as_tibble", "tibble")
import::from(tidyr, "gather")

# load the data from the gds object
wheat_data <- snpgds_parse(phys_gds)

# get the clusters
cluster <- read_rds(hdbscan)$cluster

index_chrs_chrw_csws <- which(
  (wheat_data$sample$annot$pheno == "HRS" & cluster == 5)
  | (wheat_data$sample$annot$pheno == "HRW" & cluster == 1)
  | (wheat_data$sample$annot$pheno == "SWS" & cluster == 2)
)
cpheno <- wheat_data$sample$annot$pheno[index_chrs_chrw_csws] %>% factor()
geno <- wheat_data$genotype[, index_chrs_chrw_csws]

major_allele_freq_pops <- apply(geno, 1, function(marker) {
  total_geno <- table(marker)
  major <- which.max(total_geno) %>% names()
  by(marker, cpheno, function(pop) {
    pop_geno <- table(pop)
    if (all(c("0", "2") %in% names(pop_geno))) {
      pop_geno[[major]] / sum(pop_geno[["0"]], pop_geno[["2"]])
    } else if (major %in% names(pop_geno)) {
      1
    } else {
      0
    }
  }) %>% rbind()
}) %>% t() %>% as_tibble()
colnames(major_allele_freq_pops) <- c("chrs", "chrw", "csws")

group <- "chrs_chrw_csws"
# add the genes positons to the regions table
wheat_data$snp <- wheat_data$snp %>%
  add_column(group := read_rds(josts_d)) %>%
  gather(group, D, group) %>%
  cbind(major_allele_freq_pops %>% round(4)) %>%
  rowwise() %>%
  mutate(class =
    ifelse(
      D > 0.32,
      c("CHRSD", "CHRWD", "CSWSD")[
        which.max(
          c(
            sum(abs(chrs - chrw), abs(chrs - csws)),
            sum(abs(chrw - chrs), abs(chrw - csws)), 
            sum(abs(csws - chrs), abs(csws - chrw))
          )
        )
      ],
      "None"
    )
  )
    

# overall freqs
summary(wheat_data$snp$D)
sum(wheat_data$snp$D > mean(wheat_data$snp$D)) / length(wheat_data$snp$D)

# genome freqs
for (genome in c("A", "B", "D")) {
  median(wheat_data$snp$D[which(grepl(genome, wheat_data$snp$chrom))]) %>%
    round(4) %>% str_c(genome, " median = ", .) %>% print()
}

# chromsome group freqs
for (chr_group in 1:7) {
  median(wheat_data$snp$D[which(grepl(chr_group, wheat_data$snp$chrom))]) %>%
    round(4) %>% str_c("Chr ", chr_group, " median = ", .) %>% print()
}

# group freq and median values
classes <- c("CHRSD", "CHRWD", "CSWSD", "None")

for (class in classes) {
  (sum(wheat_data$snp$class == class) / length(wheat_data$snp$class) * 100) %>%
    round(2) %>% str_c("Class ", class, " percent = ", .) %>% print()
  median(wheat_data$snp$D[which(wheat_data$snp$class == class)]) %>%
    round(2) %>% str_c("Class ", class, " median = ", .) %>% print()
}

# # check out the trend in differences between nearby Jost's D values
# iidd <- by(wheat_data$snp, wheat_data$snp$chrom, function(chrom) {
#     ret <- tibble(
#       id_a = character(), id_b = character(),
#       dists = numeric(), diffs = numeric()
#     )
#     for (i in 1:nrow(chrom)) {
#     nearby <- which(
#       chrom$pos_mb < (chrom$pos_mb[i] + 0.001) & chrom$pos_mb > chrom$pos_mb[i]
#     )
#     if (length(nearby)) {
#       ret <- ret %>% add_row(
#         id_a = chrom$id[i],
#         # a = rep(chrom$id[i], length(nearby)),
#         id_b = chrom$id[nearby],
#         dists = chrom$pos_mb[nearby] - chrom$pos_mb[i],
#         diffs = chrom$D[nearby] - chrom$D[i]
#       )
#     }
#   }
#   ret
# }) %>% do.call(rbind, .)
# iidd$diffs <- iidd$diffs %>% abs()
# nrow(iidd)
# num_markers <- unique(c(iidd$id_a, iidd$id_b)) %>% length()
# xtrm_diff <- iidd[which(iidd$diffs > 0.5), ]
# num_xtrm_markers <- unique(c(xtrm_diff$id_a, xtrm_diff$id_b)) %>% length()
# num_xtrm_markers / num_markers
# dif_freq <- sum(iidd$diffs > 0.5) / nrow(iidd)
# dif_freq * nrow(iidd)
# dists_diffs %>% ggplot(aes(dists, diffs)) +
#   geom_point()

# load the gene positions
pheno_genes <- load_genes(
  file.path(blast, "selected_pheno.csv"), base = 1
) %>% mutate(gene_type = "Phenotype Genes")
resi_genes <- load_genes(
  file.path(blast, "selected_resi.csv"), base = 1
) %>% mutate(gene_type = "Resistance Genes")

wheat_data$snp <- wheat_data$snp %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, pos_mb) %>%
  mutate(type = pmin(gene_type, class, na.rm = TRUE)) %>%
  select(-c(gene_type, class))

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
        wheat_data$max_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ] / 1e6
      ) +
      geom_point(
        aes(pos_mb, D, colour = type), shape = 16, size = 1
      ) +
      geom_smooth(
        aes(pos_mb, D), colour = colour_set[22], method = "loess", span = 0.06,
        size = 0.475, se = FALSE
      ) +
      geom_point(
        aes(pos_mb, base, colour = type), size = 1.5, shape = 25
      ) +
      scale_colour_manual(
        legend_title, values = colours_comps_genes,
        limits = levels(as.factor(wheat_data$snp$type)),
        guide = guide_legend(
          override.aes = list(
            shape = c(rep(16, 4), 25, 25)
          )
        )
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
png(
  file.path("results", "clustered_phenos_markers_josts_ds_with_genes.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()

################################################################################
# print out all markers on each chromosome
blah <- by(wheat_data$snp %>% select(-c(base, group)), wheat_data$snp$chrom,
  function(chrom) {
    write_csv(
      chrom, 
      file.path(josts_d_by_chrom,
        str_c(
          chrom$chrom[1], 
          "_all_markers_Ds_and_major_allele_freqs_by_pop_with_genes.csv"
        )
      )
    )
  }
)

################################################################################
import::from(ggrepel, "geom_text_repel")
import::from(ggplot2, "element_text", "guides", "labs")
legend_title <- "Marker Types & Genes"
# png("Results/loci/D/comps_VRN-A1.png",
#   family = "Times New Roman", width = 100, height = 62, pointsize = 1,
#   units = "mm", res = 300
# )
wheat_data$snp[
  which(
    wheat_data$snp$chrom == "2A"
    & wheat_data$snp$pos_mb > 703
    & wheat_data$snp$pos_mb < 719
  ),
] %>%
# arrange(pos_mb) %>% print(n = Inf)
  ggplot() +
    geom_point(aes(pos_mb, D, colour = type), size = 1, alpha = 0.5) +
    geom_text_repel(
      aes(pos_mb, base, colour = type, label = id),
      vjust = 1, size = 2, fontface = "bold",
      nudge_y = -0.05,
      show.legend = FALSE
    ) + 
    scale_colour_manual(
      legend_title, values = colours_comps_genes,
      limits = levels(as.factor(wheat_data$snp$type))
    ) +
    labs(x = "Postion in Mb", y = "Jost's D") +
    theme(
      legend.position = "bottom",
      # legend.text = element_text(size = 5),
      # legend.title = element_text(size = 5),
      text = element_text(size = 7.5),
      ) +
    guides(colour = guide_legend(nrow = 3, byrow = TRUE))
# dev.off()

################################################################################
# mutate(class =
#   ifelse(
#     (
#       (chrs >= 0.5 && chrw >= 0.5 && csws < 0.5)
#       || (chrs < 0.5 && chrw < 0.5 && csws >= 0.5)
#     )
#     && D > 0.29,
#     "CSWSD",
#     ifelse(
#       (
#         (chrs < 0.5 && chrw >= 0.5 && csws < 0.5)
#         || (chrs >= 0.5 && chrw < 0.5 && csws >= 0.5)
#       )
#       && D > 0.29,
#       "CHRWD",
#       ifelse(
#         (
#           (chrs >= 0.5 && chrw < 0.5 && csws < 0.5)
#           || (chrs < 0.5 && chrw >= 0.5 && csws >= 0.5)
#         )
#         && D > 0.29,
#         "CHRSD",
#         "None"
#       )
#     )
#   ),
# )
# mutate(class =
#   ifelse(
#     (chrs >= 0.6 && chrw >= 0.6 && csws <= 0.4) ||
#     (chrs <= 0.4 && chrw <= 0.4 && csws >= 0.6),
#     "CSWS vs CHRS & CHRW",
#     ifelse(
#       (chrs > 0.4 && chrs < 0.6) &&
#       ((chrw >= 0.6 && csws <= 0.4) || (chrw <= 0.4 && csws >= 0.6)),
#       "CHRW vs CSWS",
#       ifelse(
#         (chrs <= 0.4 && chrw >= 0.6 && csws <= 0.4) ||
#         (chrs >= 0.6 && chrw <= 0.4 && csws >= 0.6),
#         "CHRW vs CHRS & CSWS",
#         ifelse(
#           (chrw > 0.4 && chrw < 0.6) &&
#           ((chrs >= 0.6 && csws <= 0.4) || (chrs <= 0.4 && csws >= 0.6)),
#           "CHRS vs CSWS",
#           ifelse(
#             (chrs >= 0.6 && chrw <= 0.4 && csws <= 0.4) ||
#             (chrs <= 0.4 && chrw >= 0.6 && csws >= 0.6),
#             "CHRS vs CHRW & CSWS",
#             ifelse(
#               (csws > 0.4 && csws < 0.6) &&
#               ((chrs >= 0.6 && chrw <= 0.4) || (chrs <= 0.4 && chrw >= 0.6)),
#               "CHRS vs CHRW",
#               "None"
#             )
#           )
#         )
#       )
#     )
#   ),
#   type = pmin(gene_type, class, na.rm = TRUE)
# ) %>%
