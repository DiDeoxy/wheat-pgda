library(SNPRelate)
library(scrime)
library(poppr)
library(extrafont)

## Minimim Spanning Network stuff
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_imputed_subset.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))

## making the distance object
wheat.dist <- as.dist(1-snpgdsIBS(wheat)$ibs)
snpgdsClose(wheat)

## making the msn
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_imputed_subset_genlight.RData")
msn <- poppr.msn(wheat.genlight, distmat = wheat.dist, showplot = F) 

## plotting it
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\msn\\HRA_msn_dend.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
set.seed(106)
plot_poppr_msn(wheat.genlight, msn, inds = "which(wheat.genlight$pop ==)",
               palette = colors()[c(554, 53, 25, 26)])
title("Fig 3: Minimum Spanning Network\nof Wheat Varieties", out = T, cex = 1.5, line = -2)
dev.off()