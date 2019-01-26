library(SNPRelate)
library(vegan)
library(extrafont)
library(MASS)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\R\\functions\\colour_sets.R")
source("Analysis\\R\\functions\\funcs_draw_dend.R")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

## Minimim Spanning Network stuff
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
## making the distance object
dist <- as.dist(1-snpgdsIBS(wheat, autosome.only = F)$ibs)
snpgdsClose(wheat)

## loading the group data
load("Data\\Intermediate\\dbscan\\wheatHdbscan.RData")
# wheatHdbscan$cluster[which(wheatHdbscan$cluster == 0)] <- wheatHdbscan$cluster[which(wheatHdbscan$cluster == 0)] + 5
bests <- factor(wheatHdbscan$cluster + 1)
levels(bests) <- c("Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

## making the msn
msn <- spantree(dist)
msn$labels <- sample.id
coph <- cophenetic(msn)
conf <- scores(sammon(coph, magic = 0.35))

## plotting it
png("Results\\msn\\full_msn_hdbscan.png",
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 500)
par(oma = c(1, 1, 1, 1), mar = c(0, 0, 0, 0))
plot(msn, col = coloursDBSCAN[bests], ylab = "", xlab = "", xaxt = "n", yaxt = "n", cex = 0.4, font = 2, type = "t", ord = conf)
legends("bottomright", bests, coloursDBSCAN, "ICL Groups", cex = 0.75)
title("Minimum Spanning Network of Wheat Varieties\nColoured by ICL Group", 
      out = F, cex = 1.5, line = -3)
dev.off()

png("Results\\msn\\full_msn_desig.png",
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 500)
par(oma = c(1, 1, 1, 1), mar = c(0, 0, 0, 0))
plot(msn, col = coloursDesig[desig], ylab = "", xlab = "", xaxt = "n", yaxt = "n", cex = 0.4, font = 2, type = "t", ord = conf)
legends("bottomright", desig, coloursDesig, "Designations", cex = 0.75)
title("Minimum Spanning Network of Wheat Varieties\nColoured by Designation", 
      out = T, cex = 1.5, line = -3)
dev.off()
