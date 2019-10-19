library(SNPRelate)
library(extrafont)
library(ape)
library(dendextend)
library(RColorBrewer)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample.id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))

## PCA
pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
pc.percent <- pca$varprop*100

## lables
neighbour <- nj(1 - snpgdsIBS(wheat, autosome.only = F)$ibs)
snpgdsClose(wheat)
labels(neighbour) <- sample.id
neighbour <- root(neighbour, which(sample.id == "Selkirk"), r = TRUE)
neighbour <- chronos(neighbour, lambda = 0, model = "correlated")
neighbour <- as.dendrogram(neighbour)
neighbour <- color_branches(neighbour, k = 5, col =  c(brewer.pal(4, "Dark2"), "black"))
pop.code <- as.factor(get_leaves_branches_col(neighbour))[match(pca$sample.id, labels(neighbour))]
colours <- levels(pop.code)


## Plotting the first four dimenions against each other
png("Results\\pca\\full_pca_NJ.png", family="Times New Roman",
    width = 6, height = 6.5, pointsize = 10, units = "in", res = 300)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col = colours[as.numeric(pop.code)],
      labels=lbls, pch = 19, cex = 0.8, 
      # main = "Fig 3: First Three PCs Plotted Against Each Other\nWith Varieties Coloured by 6 K-means Groups"
      oma = c(6,2,2,2)
      )
par(xpd = T)
legend("bottom", legend=c("Group A", "Group B", "Group C", "Group D"),
       pch=19, col=brewer.pal(4, "Dark2"), cex = 0.6, horiz = T)
dev.off()

## 3D plot
# plot3d(pca$eigenvect[,1:3], col = as.integer(pop.code), type = "s", size = 1, xlab = "Eigenvector 1",
#        ylab = "Eigenvector 2", zlab = "Eigenvector 3")