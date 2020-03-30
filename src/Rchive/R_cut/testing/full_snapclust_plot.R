library(SNPRelate)
library(extrafont)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))

## PCA
pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
snpgdsClose(wheat)
pc.percent <- pca$varprop*100

load("Data\\Intermediate\\full_snapclust.Rdata")
pop.code <- sc.wheat$group

## colour pallete for pop groups
colours <- colors()[c(554, 100, 114, 53, 654, 115)]
# pie(rep(1, 6), col = colours)
palette(colours)

## Plotting the first four dimenions against each other
png("Results\\pca\\full_pca_sc.png", family="Times New Roman",
    width = 6, height = 6.5, pointsize = 10, units = "in", res = 300)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:3], col = colours[pop.code],
      labels=lbls, pch = 19, cex = 0.8, 
      # main = "Fig 3: First Three PCs Plotted Against Each Other\nWith Varieties Coloured by 6 K-means Groups"
      oma = c(6,2,2,2)
)
par(xpd = T)
legend("bottom", legend=c("Group 1", "Group 2", "Group 3",
                          "Group 4", "Group 5", "Group 6"),
       pch=19, col=1:6, cex = 0.6, horiz = T)
dev.off()

## 3D plot
# plot3d(pca$eigenvect[,1:3], col = as.integer(pop.code), type = "s", size = 1, xlab = "Eigenvector 1",
#        ylab = "Eigenvector 2", zlab = "Eigenvector 3")