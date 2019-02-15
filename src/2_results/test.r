library(plyr)
library(tidyverse)
library(RColorBrewer)
# install.packages("reshape2")
library(reshape2)
library(GGally)
library(ggrepel)


# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/colour_sets.R")
source("src/R_functions/funcs_locus_by_locus.R")

# # find the max position of any marker on each genome for xlims
# max_genome_lengths <- calc_max_genome_lengths(wheat_data)

# wheat_data <- parse_gds("mr_pruned_phys_sample_subset")
# cluster <- read_rds("Data/Intermediate/hdbscan/wheat_hdbscan.rds")$cluster

# indices <- list(
#   chrs = which(wheat_data$sample$annot$pheno == "HRS" & cluster == 5),
#   chrw = which(wheat_data$sample$annot$pheno == "HRW" & cluster == 1),
#   csws = which(wheat_data$sample$annot$pheno == "SWS" & cluster == 2)
# )
# indices[["chrs_chrw_csws"]] <- indices %>% unlist() %>% unique()

# wheat_gds <- snpgdsOpen("Data/Intermediate/GDS/mr_pruned_phys_sample_subset.gds")
# blah <- by(wheat_data$snp, wheat_data$snp$chrom, function (chrom) {
#   snp_index <- which(chrom$pos_mb > 567 & chrom$pos_mb < 607)
#   par(mfrow = c(2, 2))
#   gap_pos <- vector()
#   plots <- list()
#   if (chrom$chrom[1] == "5A") {
#     prev <- chrom$pos_mb[snp_index][1]
#     dist <- 2
#     for (cur in chrom$pos_mb[snp_index]) {
#       if (cur - prev > dist) {
#         gap_pos <- c(
#           gap_pos,
#           # finding the average positions between two markers that are
#           # ~ dist apart
#           cbind(seq(prev, cur, dist), rev(seq(cur, prev, -dist))) %>% rowMeans()
#         )
#       }
#       prev <- cur
#     }
#     print(chrom$id[snp_index])
#     for (group in names(indices)) {    
#       mat <- snpgdsLDMat(
#         wheat_gds, method = "composite", snp.id = chrom$id[snp_index],
#         sample.id = wheat_data$sample$id[indices[[group]]],
#         slide = -1
#       )$LD %>% abs()
#       # mat[is.nan(mat)] <- -.1
#       long_mat <- mat %>% melt()
#       plots[[group]] <- ggplot(long_mat, aes(Var1, Var2, fill = value)) +
#         geom_raster() +
#         scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(5, "RdYlBu")))(100))
#       # col_gaps_mat <- mat %>%
#       #   cbind(pos_mb = chrom$pos_mb[snp_index]) %>% 
#       #   as.tibble() %>%
#       #   full_join(tibble(pos_mb = gap_pos)) %>%
#       #   arrange(pos_mb) %>%
#       #   select(-pos_mb)
#       # both_gaps_mat <- t(col_gaps_mat) %>%
#       #   cbind(pos_mb = chrom$pos_mb[snp_index]) %>%
#       #   as.tibble() %>%
#       #   full_join(tibble(pos_mb = gap_pos)) %>%
#       #   arrange(pos_mb)
#       # image(
#       #   list(
#       #     x = both_gaps_mat$pos_mb, y = both_gaps_mat$pos_mb,
#       #     z = both_gaps_mat %>% select(-pos_mb) %>% as.matrix()
#       #   ),
#       #   col = colorRampPalette(rev(brewer.pal(5, "RdYlBu")))(100)
#       # )
#     }
#   }
#   plots
# })
# snpgdsClose(wheat_gds)

# plots_matrix <- ggmatrix(
#   blah[["5A"]], nrow = 2, ncol = 2,
#   legend = c(1, 1)
# )

# plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")

################################################################################
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
    ifelse(
      (chrs >= 0.6 && chrw >= 0.6 && csws <= 0.4) ||
      (chrs <= 0.4 && chrw <= 0.4 && csws >= 0.6),
      "CSWS vs CHRS & CHRW",
      ifelse(
        (chrs > 0.4 && chrs < 0.6) &&
        ((chrw >= 0.6 && csws <= 0.4) || (chrw <= 0.4 && csws >= 0.6)),
        "CHRW vs CSWS",
        ifelse(
          (chrs <= 0.4 && chrw >= 0.6 && csws <= 0.4) ||
          (chrs >= 0.6 && chrw <= 0.4 && csws >= 0.6),
          "CHRW vs CHRS & CSWS",
          ifelse(
            (chrw > 0.4 && chrw < 0.6) &&
            ((chrs >= 0.6 && csws <= 0.4) || (chrs <= 0.4 && csws >= 0.6)),
            "CHRS vs CSWS",
            ifelse(
              (chrs >= 0.6 && chrw <= 0.4 && csws <= 0.4) ||
              (chrs <= 0.4 && chrw >= 0.6 && csws >= 0.6),
              "CHRS vs CHRW & CSWS",
              ifelse(
                (csws > 0.4 && csws < 0.6) &&
                ((chrs >= 0.6 && chrw <= 0.4) || (chrs <= 0.4 && chrw >= 0.6)),
                "CHRS vs CHRW",
                "None"
              )
            )
          )
        )
      )
    ),
    type = pmin(gene_type, comp_type, na.rm = TRUE)
  ) %>%
  # mutate(comp_type = 
  #   ifelse(
  #     (chrs >= 0.5 && chrw >= 0.5 && csws < 0.5) ||
  #     (chrs < 0.5 && chrw < 0.5 && csws >= 0.5),
  #     "CSWS vs CHRS & CHRW",
  #     ifelse(
  #       (chrs < 0.5 && chrw >= 0.5 && csws < 0.5) ||
  #       (chrs >= 0.5 && chrw < 0.5 && csws >= 0.5),
  #       "CHRW vs CHRS & CSWS",
  #       ifelse(
  #         (chrs >= 0.5 && chrw < 0.5 && csws < 0.5) ||
  #         (chrs < 0.5 && chrw >= 0.5 && csws >= 0.5),
  #         "CHRS vs CHRW & CSWS",
  #         "None"
  #       )
  #     )
  #   ),
  #   type = pmin(gene_type, comp_type, na.rm = TRUE)
  # ) %>%
  select(-c(gene_type, comp_type))

legend_title <- "Marker Types & Genes"
png("Results/loci/D/comps_VRN-A1.png",
  family = "Times New Roman", width = 100, height = 62, pointsize = 1,
  units = "mm", res = 300
)
wheat_data$snp[
  which(
    wheat_data$snp$chrom == "1A"
    & wheat_data$snp$pos_mb > 8.2
    & wheat_data$snp$pos_mb < 11.6
  ),
] %>%
# arrange(pos_mb) %>% print(n = Inf)
  ggplot() +
    geom_point(aes(pos_mb, D, colour = type), size = 1) +
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
dev.off()
