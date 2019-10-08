library(SNPRelate)
library(scrime)
library(poppr)
library(extrafont)

## Minimim Spanning Network stuff
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))

## making the distance object
wheat.dist <- as.dist(1-snpgdsIBS(wheat)$ibs)
snpgdsClose(wheat)

## making the msn
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_subset_both_genlight.RData")
msn <- poppr.msn(wheat.genlight, distmat = wheat.dist, showplot = F) 

## plotting it
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\msn\\full_msn.png",
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
set.seed(111)
plot_poppr_msn(wheat.genlight, msn, inds = "", palette = colors()[c(554, 100, 114, 652, 53, 115)])
title("Fig 3: Minimum Spanning Network\nof Wheat Varieties", out = T, cex = 1.5, line = -2)
dev.off()
