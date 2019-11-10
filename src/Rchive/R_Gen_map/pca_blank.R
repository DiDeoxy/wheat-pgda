library(SNPRelate)
library(extrafont)

## loading the data
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset.gds")

## PCA
pca <- snpgdsPCA(wheat.subset, num.thread = 4, autosome.only = F)
snpgdsClose(wheat.subset)
pc.percent <- pca$varprop*100

## Plotting the first four dimenions against each other
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\pca\\pca2_blank_1_to_4.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], labels=lbls, pch = 19, cex = 0.9,
      main = paste("PCA of Pruned SNP and Sample Set"))
dev.off()