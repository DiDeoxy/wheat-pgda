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

load("Data\\Intermediate\\dbscan\\wheatHdbscan.RData")
bests <- factor(wheatHdbscan$cluster + 1)

## Plotting the first four dimenions against each other
panelfun <- function(x, y, cluster = bests, prob = wheatHdbscan$membership_prob, ...) {
  points(x, y,
         col = coloursDBSCAN[cluster],
         pch = ifelse(cluster == 1, 4, 1), # Mark noise as star
         cex = ifelse(cluster == 1, 0.5, 0.75))
  colours <- sapply(1:length(cluster), 
                   function(i) { adjustcolor(coloursDBSCAN[(cluster)[i]], alpha.f = prob[i]) })
  points(x, y, col=colours, pch=20)
}
png("Results\\pca\\full_pca_dbscan.png", family = "Times New Roman",
    width = 6, height = 7, pointsize = 10, units = "in", res = 300)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits = 2), "%", sep="")
pairs(pca$eigenvect[,1:3], 
      panel = panelfun,
      labels = lbls,
      main = "First Three PCs Plotted Against Each Other\nWith Varieties Coloured by HDBSCAN Clustering",
      oma = c(6,2,8,2)
)
par(xpd = T)
legend("bottom", legend = c("Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"),
       pch = 19, col = coloursDBSCAN, cex = 0.6, horiz = T)
dev.off()

# # 3D plot
# plot3d(pca$eigenvect[,1:3], col = coloursDBSCAN[wheatHdbscan$cluster + 1], type = "s", size = 1, xlab = "Eigenvector 1",
#        ylab = "Eigenvector 2", zlab = "Eigenvector 3")