library(SNPRelate)
library(tidyverse)
library(magrittr)
library(GGally)

# identify and remove from the overall dataset LD pruned markers
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")

# performs LD pruning to produce a set of SNPs that maximally represent the
# diversity of the genome with as little redundant info as possible used for
# estimation of relationships genome wide (i.e. clustering)
# applies multiple values of maf, bps, and ld, visual inspection of the
# resulting PCA plots identified the best subset
full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_subset_sample.gds")
for (maf in seq(0.00, 0.25, 0.05)) {
  for (bps in c(1e5, 1e6, 1e7, 1e8)) {
    for (ld in seq(0.3, 0.8, 0.1)) {
      set.seed(1000)
      snpset_ids_list <- snpgdsLDpruning(
        full,
        autosome.only = F, missing.rate = 0.1, maf = maf,
        ld.threshold = ld, slide.max.bp = bps
      )
      snps_kept <- snp_id[match(unlist(snpset_ids_list), snp_id)]

      # PCA
      pca <- snpgdsPCA(full, snp.id = snps_kept, autosome.only = F)

      png(str_c(
        "Results\\pca\\tests\\pca_LD_", ld * 10, "_bps_", bps, "_maf_",
        maf * 100, "_", length(snps_kept), ".png"
      ),
      family = "Times New Roman", width = 200, height = 200, pointsize = 5,
      units = "mm", res = 300
      )
      pairs(pca$eigenvect[, 1:3], main = str_c(
        "pca_LD_", ld * 10, "_bps_",
        bps, "_maf_", maf * 100, "_", length(snps_kept)
      ))
      dev.off()

      # dend
      full_dist <- as.dist(1 - snpgdsIBS(full, snp.id = snps_kept,
        autosome.only = F)$ibs)

      upgma_dend <- as.dendrogram(hclust(full_dist), method = "average")
      upgma_dend <- color_branches(upgma_dend, k = 5, col = colours_dend)
      label_order <- order.dendrogram(upgma_dend)

      png(str_c(
        "Results\\dend\\tests\\dend_LD_", ld * 10, "_bps_", bps, "_maf_",
        maf * 100, "_", length(snps_kept), ".png"
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
                    sample_id[label_order], facing = "clockwise",
                    niceFacing = T, cex = 0.25, adj = c(0, -0.2), font = 2)
      }, track.height = 0.15, bg.border = NA)

      draw_rects(mc, colours_mc, label_order, colors()[105])
      draw_rects(desig, colours_desig, label_order, colors()[525])

      max_height <- max(attr(upgma_dend, "height"))
      circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
        circos.dendrogram(upgma_dend, max_height = max_height)
      }, track.height = 0.3, bg.border = NA)

      circos.clear()

      title(main = paste(str_c("dend_LD_", ld * 10, "_bps_", bps, "_maf_",
        maf * 100, "_", length(snps_kept))),
            cex.main = 0.7)

      dev.off()
    }
  }
}
snpgdsClose(full)