library(SNPRelate)
library(extrafont)
library(dendextend)
library(rgl)
library(plyr) #using revalue function
library(devEMF)

## loading the gds of the data and pullling some attributes out
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset.gds")

## PCA
pca <- snpgdsPCA(wheat.subset, num.thread = 4, autosome.only = F)
snpgdsClose(wheat.subset)
pc.percent <- pca$varprop*100

## pop.code from dendrogram clustering
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
dend <- genotypes.pvclust %>% as.dendrogram
pop.code <- vector()
pop.code[labels(dend)] <- c(rep("Hard Red\nSpring", 238), rep("Chinese", 10), rep("Australian", 4),
                            rep("Prairie Spring\nRed", 28), rep("White Spring", 47), 
                            rep("Winter - Out", 3), rep("Winter", 58))
pop.code <- factor(pop.code)

## colour pallete for pop groups
colours <- colors()[c(258, 144, 554, 53, 114, 115, 450)]
palette(colours)

## Plotting the first four dimenions against each other
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\pca\\pca_dend.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.integer(pop.code), labels=lbls, pch = 19, cex = 0.8, 
      main = "Fig 3: First Four PCs Plotted Against Each Other\nWith Varieties Coloured by 7 Dendrogram Groups",
      oma = c(8,2,7,13), xpd = F)
legend(0.83, 0.64, legend=levels(pop.code), pch=19, col=1:nlevels(pop.code), cex = 0.7)
dev.off()

## 3D plot
plot3d(pca$eigenvect[,1:3], col = as.integer(pop.code), type = "s", size = 1, xlab = "Eigenvector 1", 
       ylab = "Eigenvector 2", zlab = "Eigenvector 3")