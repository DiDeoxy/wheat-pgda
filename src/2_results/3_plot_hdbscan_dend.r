library(tidyverse)
library(SNPRelate)
library(circlize)
library(dendextend)
library(extrafont)
library(ape)

source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/funcs_draw_dend.R")
source("src/R_functions/colour_sets.R")

## setting up the data
wheat_data <- parse_gds("phys_subset_sample_ld_pruned")

clusters <- factor(
  read_rds("Data/Intermediate/hdbscan/wheat_hdbscan.rds")$cluster)
levels(clusters) <- c("Noise", "Cluster 1", "Cluster 2", "Cluster 3",
                      "Cluster 4", "Cluster 5")

wheat_gds <- snpgdsOpen(
  "Data/Intermediate/GDS/phys_subset_sample_ld_pruned.gds")
## making the distance object
ibs_dist <- as.dist(1 - snpgdsIBS(wheat_gds, autosome.only = F)$ibs)
snpgdsClose(wheat_gds)

upgma_dend <- as.dendrogram(hclust(ibs_dist), method = "average")
upgma_dend <- color_branches(upgma_dend, k = 7, col = colours_dend)
label_order <- order.dendrogram(upgma_dend)

## drawing the circos plot
png("Results/dend/UPGMA_hdbscan_annot.png",
    family = "Times New Roman", width = 420, height = 420, pointsize = 30,
    units = "mm", res = 500
    )
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5,
           track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, length(label_order)), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:length(label_order), rep(0, length(label_order)),
    wheat_data$sample$id[label_order], facing = "clockwise", niceFacing = T,
    cex = 0.25, adj = c(0, -0.2), font = 2)
},
track.height = 0.15, bg.border = NA)

draw_rects(
  wheat_data$sample$annot$era, colours_era, label_order, colors()[468])
draw_rects(
  wheat_data$sample$annot$bp, colours_bp, label_order, colors()[383])
draw_rects(
  wheat_data$sample$annot$mc, colours_mc, label_order, colors()[105])
draw_rects(
  wheat_data$sample$annot$pheno, colours_pheno, label_order, colors()[525])
draw_rects(
  clusters, colours_hdbscan, label_order, colors()[109])

max_height <- max(attr(upgma_dend, "height"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(upgma_dend, max_height = max_height)
}, track.height = 0.3, bg.border = NA)

circos.clear()

## Legends
pch <- 19
cex <- 0.45
legend("topright", legend = levels(wheat_data$sample$annot$era),
  title = "Period of Release", pch = pch, col = colours_era, cex = cex,
  bg = colors()[468])
legend("bottomright", legend = levels(wheat_data$sample$annot$bp),
  title = "Breeding Program/Origin", pch = pch, col = colours_bp, cex = 0.35,
  bg = colors()[383])
legend("bottomleft", legend = levels(wheat_data$sample$annot$mc),
  title = "Market Class", pch = pch, col = colours_mc, cex = cex,
  bg = colors()[105])
legend("topleft", legend = levels(wheat_data$sample$annot$pheno),
  title = "Phenotype", pch = pch, col = colours_pheno, cex = cex,
  bg = colors()[525])
legend("center", legend = levels(clusters), title = "HDBSCAN Clusters",
  pch = pch, col = colours_hdbscan_legend, cex = cex, bg = colors()[109])

title(
  main = str_c(
    "UPGMA Dendrogram of 364 Varieties\nwith HDBSCAN Clusters and ",
    "Categorical Data In Surrounding Rows"
  ),
  cex.main = 0.7
)
dev.off()