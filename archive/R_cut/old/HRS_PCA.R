library(SNPRelate)
library(extrafont)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("Data\\Formatted\\HRS_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))

## PCA
pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
snpgdsClose(wheat)
pc.percent <- pca$varprop*100

load("Data\\Formatted\\hrs_kmeans.RDATA")
pop.code <- bests

## colour pallete for pop groups
colours <- colors()[c(554, 100, 654, 53)]
palette(colours)
# pie(rep(1,4), col = colours)

## Plotting the first four dimenions against each other
png("Analysis\\Figures\\pca\\hrs_pca.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
# lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
# pairs(pca$eigenvect[,1:3], col = colours[pop.code],
#       labels=lbls, pch = 19, cex = 0.8, 
#       main = "Fig 3: First Three PCs Plotted Against Each Other\nWith Varieties Coloured by 6 K-means Groups"
#       )
plot(pca$eigenvect[,1:2], col = colours[pop.code], pch = 19, cex = 0.8, 
      # main = "Fig 3: First Three PCs Plotted Against Each Other\nWith Varieties Coloured by 6 K-means Groups"
)
# legend(0.8, 0.65, legend=unique(pop.code), pch=19, col=1:nlevels(pop.code), cex = 0.6)
dev.off()

## 3D plot
# plot3d(pca$eigenvect[,1:3], col = as.integer(pop.code), type = "s", size = 1, xlab = "Eigenvector 1",
#        ylab = "Eigenvector 2", zlab = "Eigenvector 3")