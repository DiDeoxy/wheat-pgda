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
  arrange(chrom, pos_mb) %>%
  rowwise() %>%
  # mutate(comp_type = 
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
  #   type = pmin(gene_type, comp_type, na.rm = TRUE)
  # ) %>%
  mutate(comp_type = 
    ifelse(
      (chrs >= 0.5 && chrw >= 0.5 && csws < 0.5) ||
      (chrs < 0.5 && chrw < 0.5 && csws >= 0.5),
      "CSWS vs CHRS & CHRW",
      ifelse(
        (chrs < 0.5 && chrw >= 0.5 && csws < 0.5) ||
        (chrs >= 0.5 && chrw < 0.5 && csws >= 0.5),
        "CHRW vs CHRS & CSWS",
        ifelse(
          (chrs >= 0.5 && chrw < 0.5 && csws < 0.5) ||
          (chrs < 0.5 && chrw >= 0.5 && csws >= 0.5),
          "CHRS vs CHRW & CSWS",
          "None"
        )
      )
    ),
    type = pmin(gene_type, comp_type, na.rm = TRUE)
  ) %>%
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
        aes(pos_mb, D, colour = type), shape = 16, size = 1
      ) +
      geom_smooth(
        aes(pos_mb, D), colour = colour_set[22], method = "loess", span = 0.06,
        size = 0.475, se = FALSE
      ) +
      geom_point(
        aes(pos_mb, base, colour = type), size = 1.5, shape = 25
      ) +
      # geom_text_repel(
      #   aes(pos_mb, base, colour = type, label = id), angle = 90, hjust = 0,
      #   vjust = -1, size = 3, fontface = "bold", nudge_y = -0.07,
      #   nudge_x =
      #     ifelse(chrom %in% c("1A", "1D"), 80,
      #       ifelse(chrom %in% c("2D", "4A"), -60, 40)
      #     ),
      #   show.legend = FALSE
      # ) + 
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
png("Results/loci/D/comps_D_test.png",
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()

################################################################################
# print out all markers on each chromosome
base <- "Results/loci/D/by_chrom"
ifelse(! dir.exists(file.path(base)), dir.create(file.path(base)), FALSE)
blah <- by(wheat_data$snp %>% select(-c(base, group)), wheat_data$snp$chrom, function (chrom) {
  write_csv(chrom, file.path(base, str_c(chrom$chrom[1], "_all_markers_Ds_and_major_allele_freqs_by_pop_with_genes.csv")))
})

################################################################################
# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

wheat_data$snp <- wheat_data$snp %>%
  add_group_stat(group) %>%
  gather(group, D, group) %>%
  cbind(major_allele_freq_pops %>% round(4)) %>%
  rowwise() %>%
  mutate(type = 
    ifelse(
      (chrs >= 0.5 && chrw >= 0.5 && csws < 0.5) ||
      (chrs < 0.5 && chrw < 0.5 && csws >= 0.5),
      "CSWS vs CHRS & CHRW",
      ifelse(
        (chrs < 0.5 && chrw >= 0.5 && csws < 0.5) ||
        (chrs >= 0.5 && chrw < 0.5 && csws >= 0.5),
        "CHRW vs CHRS & CSWS",
        ifelse(
          (chrs >= 0.5 && chrw < 0.5 && csws < 0.5) ||
          (chrs < 0.5 && chrw >= 0.5 && csws >= 0.5),
          "CHRS vs CHRW & CSWS",
          "None"
        )
      )
    )
  )

chrom_D_sum <- by(wheat_data$snp, wheat_data$snp$chrom, function (chrom) {
  summary(chrom$D[which(chrom$D > quantile(wheat_data$snp$D, 0.25))], na.rm = TRUE)
}) %>% do.call(rbind, .) %>% as_tibble()

sum(wheat_data$snp$D > mean(wheat_data$snp$D)) / length(wheat_data$snp$D)
mean(wheat_data$snp$D[which(wheat_data$snp$D > quantile(wheat_data$snp$D, 0.25))])
summary(wheat_data$snp$D)
sum(wheat_data$snp$D < quantile(wheat_data$snp$D, 0.25)) / nrow(wheat_data$snp)

types <- c(
  "CHRS vs CHRW & CSWS", "CHRW vs CHRS & CSWS", "CSWS vs CHRS & CHRW", "None"
)

for (type in types) {
  (sum(wheat_data$snp$type == type) / length(wheat_data$snp$D)) %>% print()
  median(wheat_data$snp$D[
    which(wheat_data$snp$type == type)
  ]) %>% print()
}



for (genome in c("A", "B", "D")) {
  # mean(wheat_data$snp$D[
  #   which(
  #     wheat_data$snp$D > quantile(wheat_data$snp$D, 0.25)
  #     & grepl(genome, wheat_data$snp$chrom)
  #   )
  # ]) %>% print()
  median(wheat_data$snp$D[which(grepl(genome, wheat_data$snp$chrom))]) %>% print()
}

for (chrom_set in 1:7) {
  median(wheat_data$snp$D[which(grepl(chrom_set, wheat_data$snp$chrom))]) %>% print()
}

jost_D_variation <- function(chrom, n) {
  vars <- vector()
  for (i in 1:(nrow(chrom) - n)) {
    vars <- c(vars, sd(chrom$D[i:(i + n)]))
  }
  vars
}


wheat_gds <- snpgdsOpen("Data/Intermediate/GDS/mr_pruned_phys_sample_subset.gds")
lds <- by(wheat_data$snp, wheat_data$snp$chrom, function (chrom) {
  ld_mat <- snpgdsLDMat(
    wheat_gds, method = "composite", snp.id = chrom$id, slide = -1
  )$LD %>% abs()
  lds <- list()
  for (i in 1:nrow(ld_mat)) {
    nearby <- which(chrom$pos_mb > (chrom$pos_mb[i] - 2) & chrom$pos_mb < (chrom$pos_mb[i] + 2))
    if (length(nearby)) {
      lds[[i]] <- cbind(
        abs(chrom$pos_mb[i] - chrom$pos_mb[nearby]), ld_mat[i, nearby] 
      )
    }
  }
  lds %>% do.call(rbind, .)
}) %>% do.call(rbind, .)
snpgdsClose(wheat_gds)
lds %>% cor(., use = "complete")

var_mn_ds <- by(wheat_data$snp, wheat_data$snp$chrom, function (chrom) {
  stats <- tibble(var = vector(), mn = vector())
  for (i in 1:nrow(chrom)) {
    nearby <- which(chrom$pos_mb > (chrom$pos_mb[i] - 20) & chrom$pos_mb < (chrom$pos_mb[i] + 20))
    stats <- stats %>% add_row(
      var = sd(chrom$D[nearby]), mn = mean(chrom$D[nearby])
    )
  }
  stats
}) %>% do.call(rbind, .)
mean(var_mn_ds$var, na.rm = TRUE) /
mean(var_mn_ds$mn, na.rm = TRUE)

with(wheat_data$snp, which(chrom == "5A" & pos_mb > 400 & pos_mb < 500)) %>% length()
with(wheat_data$snp, which(chrom == "5A" & pos_mb > 400 & pos_mb < 500 & grepl("None", type))) %>% length()

which(wheat_data$snp$chrom = "5A" & wheat)

mean(var_ds, na.rm = TRUE) / mean(wheat_data$snp$D)

mean(lds %>% unlist(), na.rm = TRUE)

################################################################################
  # mutate(comp_type = 
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
  #   )
  # ) %>%