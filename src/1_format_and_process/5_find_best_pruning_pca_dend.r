library(SNPRelate)
library(tidyverse)
library(circlize)
library(dendextend)
library(GGally)
library(extrafont)

source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/funcs_draw_dend.R")
source("src/R_functions/colour_sets.R")

wheat_data <- parse_gds("phys_subset_sample")

# performs LD pruning to produce a set of SNPs that maximally represent the
# diversity of the genome with as little redundant info as possible used for
# estimation of relationships genome wide (i.e. clustering)
# applies multiple values of maf, bps, and ld, visual inspection of the
# resulting PCA plots identified the best subset
wheat_gds <- snpgdsOpen("Data/Intermediate/GDS/phys_subset_sample.gds")
maf <- 0.05
for (max_dist in c(5e6, 1e7, 2e7, 3e7)) {
  for (ld in seq(0.8, 0.8, 0.1)) {
    # ld pruned set of markes
    set.seed(1000)
    kept_id <- unlist(snpgdsLDpruning(wheat_gds, autosome.only = FALSE,
      maf = maf, slide.max.bp = max_dist, ld.threshold = ld))

    # pca
    pca <- snpgdsPCA(wheat_gds, snp.id = kept_id, autosome.only = F)

    # pca pairs
    ggpairs_pca <- pca$eigenvect[, 1:3] %>%
      as.tibble() %>%
      rename(PC1 = V1, PC2 = V2, PC3 = V3) %>%
      ggpairs(
        title = str_c(
          "max_dist_", max_dist / 1e6, "_LD_", ld * 10, "_num_snps_",
          length(kept_id)
      )
    )

    # plot pairs
    png(
      str_c(
        "Results/pca/tests/pca_", "max_dist_", max_dist / 1e6,
        "_LD_", ld * 10, "_num_snps_", length(kept_id), ".png"
      ),
      family = "Times New Roman", width = 200, height = 200, pointsize = 12,
      units = "mm", res = 300
    )
    print(ggpairs_pca)
    dev.off()

    # dend
    full_dist <- as.dist(1 - snpgdsIBS(wheat_gds, snp.id = kept_id,
      autosome.only = F)$ibs)
    upgma_dend <- as.dendrogram(hclust(full_dist), method = "average") %>%
      color_branches(k = 7, col = colours_dend)
    label_order <- order.dendrogram(upgma_dend)

    # plot dend
    png(
      str_c(
        "Results/dend/tests/dend_", "max_dist_", max_dist / 1e6,
        "_LD_", ld * 10, "_num_snps_", length(kept_id), ".png"
      ),
      family = "Times New Roman", width = 200, height = 200, pointsize = 15,
      units = "mm", res = 300
    )
    circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5,
      track.margin = c(0.005, 0.005))
    circos.initialize("foo", xlim = c(0, length(label_order)),
      sector.width = 1)
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      circos.text(1:length(label_order), rep(0, length(label_order)),
                  wheat_data$sample$id[label_order], facing = "clockwise",
                  niceFacing = T, cex = 0.25, adj = c(0, -0.2), font = 2)
    },
      track.height = 0.15, bg.border = NA)

    draw_rects(
      wheat_data$sample$annot$mc, colours_mc, label_order, colors()[105])
    draw_rects(
      wheat_data$sample$annot$pheno, colours_pheno, label_order, colors()[525])

    max_height <- max(attr(upgma_dend, "height"))
    circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
      circos.dendrogram(upgma_dend, max_height = max_height)
    },
      track.height = 0.3, bg.border = NA)
    circos.clear()

    title(
      main = str_c(
        "dend_", "max_dist_", max_dist / 1e6, "_LD_", ld * 10, "_num_snps_",
        length(kept_id)
      ),
      cex.main = 0.7)
    dev.off()
  }
}
snpgdsClose(wheat_gds)