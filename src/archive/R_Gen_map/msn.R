library(SNPRelate)
library(scrime)
library(poppr)
library(extrafont)

## Minimim Spanning Network stuff
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset.gds")
genotypes <- read.gdsn(index.gdsn(wheat.subset, "genotype"))
snpgdsClose(wheat.subset)

genotypes <- replace(genotypes, genotypes == 3, NA)
genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes2 <- knncatimpute(genotypes)

wheat.dist <- as.dist(1-cor(genotypes2, use = "pairwise.complete.obs"))

load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genlight.RData")

set.seed(100)
msn <- poppr.msn(wheat.subset, distmat = wheat.dist, showplot = F) 

png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\msn\\msn_dend.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
set.seed(100)
plot_poppr_msn(wheat.subset, msn, inds = "peter", palette = colors()[c(258, 144, 554, 53, 114, 115, 450)])
title("Fig 3: Minimum Spanning Network\nof Wheat Varieties", out = T, cex = 1.5, line = -2)
dev.off()
