library(SNPRelate)
library(fpc)
library(extrafont)
library(dendextend)
library(rgl)

## loading the data
wheat.subset <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\wheat_subset.gds")
genotypes <- t(read.gdsn(index.gdsn(wheat.subset, "genotype")))

## PCA
pca <- snpgdsPCA(wheat.subset, num.thread = 4, autosome.only = F)
snpgdsClose(wheat.subset)
pc.percent <- pca$varprop*100

# #revaluing and imputing gentotypes
genotypes <- replace(genotypes, genotypes == 3, NA)
genotypes <- replace(genotypes, genotypes == 0, 1)
genotypes2 <- t(knncatimpute(genotypes))

## clustering
#creating the dist object
wheat.dist <- as.dist(1-cor(genotypes2, use = "pairwise.complete.obs"))
# computing likeliest pamk groups
pop.code <- pamk(wheat.dist, krange=2:10, diss = T)$pamobject$clustering

## colour pallete for PAM pop groups
colours <- colors()[c(554, 115)]
palette(colours)

## Plotting the first four dimenions against each other
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\pca\\pca_pam.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=as.integer(pop.code), labels=lbls, pch = 19, cex = 0.8,
      main = "First Four PCs Plotted Against Each Other
      With Varieties Coloured by PAM Groups", oma = c(8,2,7,13), xpd = F)
legend(0.84, 0.85, legend=c("Group 1", "Group 2"), pch = 19, col=1:2, cex = 0.7, ncol = 6)
dev.off()

## 3D plot
plot3d(pca$eigenvect[,1:3], col = as.integer(pop.code), type = "s", size = 1, xlab = "Eigenvector 1", 
       ylab = "Eigenvector 2", zlab = "Eigenvector 3")