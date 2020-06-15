# load needed packages and functions
source(file.path("src", "file_paths.R"))
suppressPackageStartupMessages(library(tidyverse))
library(SNPRelate)
library(circlize)
source(draw_dend)
library(dendextend)
library(extrafont)

# create needed file paths and directories
test_pruning <- file.path("results", "test_pruning")
ifelse(! dir.exists(test_pruning), dir.create(test_pruning), FALSE)
test_pruning_pca <- file.path(test_pruning, "pca")
ifelse(! dir.exists(test_pruning_pca), dir.create(test_pruning_pca), FALSE)
test_pruning_dend  <- file.path(test_pruning, "dend")
ifelse(! dir.exists(test_pruning_dend), dir.create(test_pruning_dend), FALSE)

# load the data
wheat_data <- parse_gds(file.path(phys_gds))

# performs LD pruning to produce a set of SNPs that maximally represent the
# diversity of the genome with as little redundant info as possible used for
# estimation of relationships genome wide (i.e. clustering)
# applies multiple values of maf, bps, and ld, visual inspection of the
# resulting PCA plots identified the best subset
wheat_gds <- snpgdsOpen(file.path(phys_gds))
for (max_dist in c(5e6, 1e7, 1.5e7, 2e7)) {
  for (ld in seq(0.5, 0.8, 0.05)) {

    # ld pruned set of markes
    set.seed(1000)
    kept_id <- unlist(
      snpgdsLDpruning(
        wheat_gds, autosome.only = FALSE, slide.max.bp = max_dist,
        ld.threshold = ld
      )
    )

    # a title for the output graphs
    title <- str_c(
      "max_dist_", max_dist / 1e6, "_LD_", ld * 100, "_num_snps_",
      length(kept_id), ".png"
    )

    # pca
    pca <- snpgdsPCA(wheat_gds, snp.id = kept_id, autosome.only = F)

    # pca pairs
    ggpairs_pca <- pca$eigenvect[, 1:3] %>%
      as.tibble() %>%
      dplyr::rename(PC1 = V1, PC2 = V2, PC3 = V3) %>%
      ggpairs(title = title)
    )

    # plot pairs
    png(
      file.path(test_pruning_pca, title),
      family = "Times New Roman", width = 200, height = 200, pointsize = 12,
      units = "mm", res = 300
    )
    print(ggpairs_pca)
    dev.off()

    # dend
    full_dist <- as.dist(1 - snpgdsIBS(wheat_gds, snp.id = kept_id,
      autosome.only = F)$ibs)
    upgma_dend <- hclust(full_dist) %>%
      as.dendrogram(method = "average") %>%
      color_branches(k = 9, col = colours_dend) %>%
      set("branches_lwd", 2)
    label_order <- order.dendrogram(upgma_dend)

    # plot dend
    png(
      file.path(test_pruning_dend, title),
      family = "Times New Roman", width = 200, height = 200, pointsize = 15,
      units = "mm", res = 300
    )
    circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5,
      track.margin = c(0.005, 0.005))
    circos.initialize("foo", xlim = c(0, length(label_order)), sector.width = 1)
    circos.track(ylim = c(0, 1), track.height = 0.15, bg.border = NA,
      panel.fun = function(x, y) {
        circos.text(
          1:length(label_order), rep(0, length(label_order)),
          wheat_data$sample$id[label_order], facing = "clockwise",
          niceFacing = T, cex = 0.25, adj = c(0, -0.2), font = 2
        )
      }
    )

    # rect rows
    draw_rects(
      wheat_data$sample$annot$mc, colours_mc, label_order, colors()[105]
    )
    draw_rects(
      wheat_data$sample$annot$pheno, colours_pheno, label_order, colors()[525]
    )

    ## Legends
    pch <- 19
    cex <- 0.45
    legend(
      "bottomleft", legend = levels(wheat_data$sample$annot$mc), box.lwd = 2,
      title = "Market Class", pch = pch, col = colours_mc, cex = cex,
      bg = colors()[105]
    )
    legend(
      "topleft", legend = levels(wheat_data$sample$annot$pheno), box.lwd = 2,
      title = "Phenotype", pch = pch, col = colours_pheno, cex = cex,
      bg = colors()[525]
    )

    max_height <- max(attr(upgma_dend, "height"))
    circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
      circos.dendrogram(upgma_dend, max_height = max_height)
    },
      track.height = 0.3, bg.border = NA)
    circos.clear()

    title(main = title, ex.main = 0.7)
    dev.off()
  }
}
snpgdsClose(wheat_gds)