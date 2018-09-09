library(SNPRelate)
library(circlize)
library(dendextend)
library(extrafont)
library(ape)

source("src\\R_functions\\funcs_draw_dend.R")
source("src\\R_functions\\colour_sets.R")

## setting up the data
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("src\\R_functions\\data_loading.R")

load("Data\\Intermediate\\dbscan\\wheatHdbscan.RData")
# wheatHdbscan$cluster[which(wheatHdbscan$cluster == 0)] <- wheatHdbscan$cluster[which(wheatHdbscan$cluster == 0)] + 5
bests <- factor(wheatHdbscan$cluster + 1)
levels(bests) <- c("Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
## making the distance object
dist <- as.dist(1-snpgdsIBS(wheat, autosome.only = F)$ibs)
snpgdsClose(wheat)

consensusTree <- as.dendrogram(hclust(dist), method = "average")
consensusTree <- color_branches(consensusTree, k = 6, col = coloursDend)
# clusters <- cutree(consensusTree, k = 6)
# # consensusTree <- as.dendrogram(wheatHdbscan$hc)
# # consensusTree <- raise.dendrogram(consensusTree, -12.5)
labelOrder <- order.dendrogram(consensusTree)

## drawing the circos plot
png("Results\\dend\\full_dbscan_clusters.png", width = 11, height = 11, pointsize = 20, units = "in", res = 500)
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, length(labelOrder)), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:length(labelOrder), rep(0, length(labelOrder)), sample.id[labelOrder],
              facing = "clockwise", niceFacing = T, cex = 0.25, adj = c(0, -0.2), font = 2)
}, track.height = 0.15, bg.border = NA)

drawRects(era, coloursEra, labelOrder, colors()[468])
drawRects(bp, coloursBp, labelOrder, colors()[383])
drawRects(mc, coloursMc, labelOrder, colors()[105])
drawRects(desig, coloursDesig, labelOrder, colors()[525])
drawRects(bests, coloursDBSCAN, labelOrder, colors()[109])

maxHeight <- max(attr(consensusTree, "height"))
circos.track(ylim = c(0, maxHeight), panel.fun = function(x, y) {
  circos.dendrogram(consensusTree, max_height = maxHeight)
}, track.height = 0.3, bg.border = NA)

circos.clear()

## Legends
legends("topright", era, coloursEra, "Period of Release", bg = colors()[468])
legends("bottomright", bp, coloursBp, "Breeding Program/Origin", cex = 0.35, bg = colors()[383])
legends("bottomleft", mc, coloursMc, "Market Class", bg = colors()[105])
legends("topleft", desig, coloursDesig, "Phenotype", bg = colors()[525])
legends("center", bests, coloursDBSCAN, "HDBSCAN Clusters", bg = colors()[109])

title(main = "Circularized UPGMA Dendrogram of Wheat Samples\nwith HDBSCAN Clusters and Variety Data In Surrounding Rows", cex.main = 0.8)

dev.off()

## testing out ggtree
# neighbour <- nj(dist)
# neighbour$tip.label <- sample.id
# neighbour <- groupOTU(neighbour, split(sample.id, bests))
# 
# tree <- ggtree(neighbour, aes(color = group)) + geom_tiplab2(align = TRUE, linesize = 0.5, cex = 1.5) +
#   scale_colour_manual(values = coloursDBSCAN2)

# png("Results\\dend\\full_dbscan_clusters_nj.png", width = 11, height = 11, pointsize = 20, units = "in", res = 500)
#   gheatmap(tree, data.frame(era), offset = 2, width = 0.5) # doesnt work mysteriously
# dev.off()
