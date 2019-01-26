library(SNPRelate)
library(extrafont)
library(dendextend)
library(rgl)

## loading the gds of the data and pullling some attributes out
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset.gds")
pop.code <- read.gdsn(index.gdsn(wheat.subset, "samp.annot/designation"))

## PCA
pca <- snpgdsPCA(wheat.subset, num.thread = 4, autosome.only = F)
snpgdsClose(wheat.subset)
pc.percent <- pca$varprop*100

## colour pallete for MC pop groups
colours <- colors()[c(554, 115, 153, 26, 153, 450, 114, 652, 153)]
palette(colours)

## Plotting the first four dimenions against each other
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\pca\\pca_designation.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.integer(pop.code), labels=lbls, 
      pch = c(rep(19, 2), 1, rep(19, 6))[as.integer(pop.code)],
      cex = 0.8, lwd = c(rep(1, 2), 0.50, rep(1, 6))[as.integer(pop.code)],
      main = "First Four PCs Plotted Against Each Other
      With Varieties Coloured by Designation", oma = c(8,2,7,13), xpd = F)
legend(0.84, 0.85, legend=levels(pop.code), pch = c(rep(19, 2), 1, rep(19, 6)), col=1:nlevels(pop.code), cex = 0.7)
dev.off()

## 3D plot
plot3d(pca$eigenvect[,1:3], col = as.integer(pop.code), type = "s", size = 1, xlab = "Eigenvector 1", 
       ylab = "Eigenvector 2", zlab = "Eigenvector 3")