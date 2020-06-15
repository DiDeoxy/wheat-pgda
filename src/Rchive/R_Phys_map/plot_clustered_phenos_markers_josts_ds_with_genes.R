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
