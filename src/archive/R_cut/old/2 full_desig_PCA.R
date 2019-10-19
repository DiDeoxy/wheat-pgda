library(SNPRelate)
library(extrafont)
# library(rgl)
# install.packages("rgl")

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\R\\functions\\colour_sets.R")

# gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
# source("Analysis\\R\\functions\\data_loading.R")

## PCA
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
snpgdsClose(wheat)

pc.percent <- pca$varprop*100

load("Data\\Intermediate\\dbscan\\dbscan.RData")
groups <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6")

## Plotting the first four dimenions against each other
png("Results\\pca\\full_pca_dbscan.png", family = "Times New Roman",
    width = 6, height = 7, pointsize = 10, units = "in", res = 300)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits = 2), "%", sep="")
pairs(pca$eigenvect[,1:3], col = coloursDBSCAN[wheatHdbscan$cluster + 1],
      labels=lbls, pch = 19, cex = 0.8, 
      main = "First Three PCs Plotted Against Each Other\nWith Varieties Coloured by phenotype",
      oma = c(6,2,8,2)
)
par(xpd = T)
legend("bottom", legend = groups,
       pch = 19, col = coloursDBSCAN, cex = 0.6, horiz = T)
dev.off()

# # 3D plot
plot3d(pca$eigenvect[,1:3], col = coloursDBSCAN[wheatHdbscan$cluster + 1], type = "s", size = 1, xlab = "Eigenvector 1",
       ylab = "Eigenvector 2", zlab = "Eigenvector 3")