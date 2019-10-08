library(SNPRelate)
library(circlize)
library(dendextend)
library(extrafont)
library(ape)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\R\\functions\\funcs_draw_dend.R")
source("Analysis\\R\\functions\\colour_sets.R")

## setting up the data
gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

load("Data\\Intermediate\\dbscan\\wheatHdbscan.RData")
wheatHdbscan$cluster[which(wheatHdbscan$cluster == 0)] <- wheatHdbscan$cluster[which(wheatHdbscan$cluster == 0)] + 5
bests <- factor(wheatHdbscan$cluster)
levels(bests) <- c("Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

consensusTree <- as.dendrogram(wheatHdbscan$hc)
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
}, track.height = 0.15, bg.border = NA)

circos.clear() 

## Legends
legends("topright", era, coloursEra, "Period of Release", bg = colors()[468])
legends("bottomright", bp, coloursBp, "Breeding Program/Origin", cex = 0.35, bg = colors()[383])
legends("bottomleft", mc, coloursMc, "Market Class", bg = colors()[105])
legends("topleft", desig, coloursDesig, "Phenotype", bg = colors()[525])
legends("center", bests, coloursDBSCAN, "dbscan Cluster", bg = colors()[109])

title(main = "Varieties Ordered By DBSCAN cluster, Phenotype,\nMarket Class, Era of Release, and Name", cex.main = 0.8)

dev.off()


