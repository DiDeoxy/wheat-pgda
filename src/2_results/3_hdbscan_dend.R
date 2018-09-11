library(tidyverse)
library(SNPRelate)
library(circlize)
library(dendextend)
library(extrafont)
library(ape)

source("src\\R_functions\\funcs_draw_dend.R")
source("src\\R_functions\\colour_sets.R")

## setting up the data
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds"
source("src\\R_functions\\data_loading.R")

clusters <- factor(
  read_rds("Data\\Intermediate\\dbscan\\full_hdbscan.rds")$cluster)
levels(clusters) <- c("Noise", "Cluster 1", "Cluster 2", "Cluster 3",
                      "Cluster 4", "Cluster 5")

full <- snpgdsOpen(
  "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds")
## making the distance object
full_dist <- as.dist(1-snpgdsIBS(full, autosome.only = F)$ibs)
snpgdsClose(full)

upgma_dend <- as.dendrogram(hclust(full_dist), method = "average")
upgma_dend <- color_branches(upgma_dend, k = 5, col = colours_dend)
label_order <- order.dendrogram(upgma_dend)

## drawing the circos plot
png("Results\\dend\\full_dbscan_clusters.png", width = 11, height = 11,
    pointsize = 20, units = "in", res = 500)
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, 
           track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, length(label_order)), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:length(label_order), rep(0, length(label_order)),
              sample_id[label_order], facing = "clockwise", niceFacing = T,
              cex = 0.25, adj = c(0, -0.2), font = 2)
}, track.height = 0.15, bg.border = NA)

draw_rects(era, colours_era, label_order, colors()[468])
draw_rects(bp, colours_bp, label_order, colors()[383])
draw_rects(mc, colours_mc, label_order, colors()[105])
draw_rects(desig, colours_desig, label_order, colors()[525])
draw_rects(clusters, colours_dbscan, label_order, colors()[109])

max_height <- max(attr(upgma_dend, "height"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(upgma_dend, max_height = max_height)
}
, track.height = 0.3, bg.border = NA)

circos.clear()

## Legends
pch <- 19
cex <- 0.45
legend("topright", legend = levels(era), title = "Period of Release",
       pch = pch, col = colours_era, cex = cex, bg = colors()[468])
legend("bottomright", legend = levels(bp), title = "Breeding Program/Origin",
       pch = pch, col = colours_bp, cex = 0.35, bg = colors()[383])
legend("bottomleft", legend = levels(mc), title = "Market Class", pch = pch,
       col = colours_mc, cex = cex, bg = colors()[105])
legend("topleft", legend = levels(desig), title = "Phenotype", pch = pch,
       col = colours_desig, cex = cex, bg = colors()[525])
legend("center", legend = levels(clusters), title = "HDBSCAN Clusters",
       pch = pch, col = colours_dbscan, cex = cex, bg = colors()[109])

title(main = paste("UPGMA Dendrogram of full Samples\nwith",
                   "HDBSCAN Clusters and Variety Data In Surrounding Rows"),
                   cex.main = 0.6)

dev.off()

## testing out ggtree
# neighbour <- nj(dist)
# neighbour$tip.label <- sample.id
# neighbour <- groupOTU(neighbour, split(sample.id, clusters))
# 
# tree <- ggtree(neighbour, aes(color = group)) + geom_tiplab2(align = TRUE, linesize = 0.5, cex = 1.5) +
#   scale_colour_manual(values = coloursDBSCAN2)

# png("Results\\dend\\full_dbscan_clusters_nj.png", width = 11, height = 11, pointsize = 20, units = "in", res = 500)
#   gheatmap(tree, data.frame(era), offset = 2, width = 0.5) # doesnt work mysteriously
# dev.off()
