library(SNPRelate)
library(fpc)
library(extrafont)

wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_subset_both.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))

# sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
# snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))

set.seed(1004)
# kmeans <- kmeans(t(genotypes), centers=3, nstart = 100)
# pamk <- pamk(1 - snpgdsIBS(wheat, autosome.only = F)$ibs, krange=3, diss = T)
# neighbour joining
neighbour <- nj(1 - snpgdsIBS(wheat, autosome.only = F)$ibs)
labels(neighbour) <- sample.id
net.struct <- read.csv("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\NetStruct\\1\\0_150.csv", stringsAsFactors = F)[,2:3]

## PCA of the different SNP sets and the full set
pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
pc.percent <- pca$varprop*100
snpgdsClose(wheat)

## colour pallete for oam pop groups
colours <- c("red2", "darkturquoise", "orange")
palette(colours)

# ## Plotting the first two dimenions against each other
# png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\pca\\HRS_kmeans.png",
#     family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
# lbls <- paste("PC", 1:2, ": ", format(pc.percent[1:2], digits=2), "%", sep="")
# plot(pca$eigenvect[,1:2], col=kmeans$cluster, pch = 19, cex = 0.8,
#      main = "Three Kmeans Groups",
#      xlab = lbls[1], ylab = lbls[2])
# # pairs(pca$eigenvect[,1:4], col=kmeans$cluster, labels = lbls, pch = 19, cex = 0.8, 
# #      main = "First Four PCs Plotted Against Each Other
# #      With Varieties Coloured by 3 PAM Groups")
# legend("topleft", legend=c("Hard Red Spring 4", "Hard Red Spring 5", "Prairie Spring Red"), pch=19, col=1:3, cex = 0.8, ncol = 2)
# dev.off()
# 
# ## Plotting the first four dimenions against each other
# png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\pca\\HRS_pam.png",
#     family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
# lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
# plot(pca$eigenvect[,1:2], col=pamk$pamobject$clustering, pch = 19, cex = 0.8, 
#      main = "Three PAM Groups",
#      xlab = lbls[1], ylab = lbls[2])
# legend("topleft", legend=c("Hard Red Spring 4", "Hard Red Spring 5", "Prairie Spring Red"), pch=19, col=1:3, cex = 0.8, ncol = 2)
# dev.off()

## Plotting the first four dimenions against each other
colours <- c("orange", "red2", "darkturquoise")
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\pca\\HRS_net_struct.png",
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
plot(pca$eigenvect[,1:2], col = sapply(net.struct[,2][match(labels(neighbour), net.struct[,1])], function(x) { if (is.na(x)) {return("grey0")} else { return(colours[x + 1]) } } ),
     pch = 19, cex = 0.8, main = "Three PAM Groups", xlab = lbls[1], ylab = lbls[2])
legend("topleft", legend=c("Hard Red Spring 4", "Hard Red Spring 5", "Prairie Spring Red"), pch=19, col=1:3, cex = 0.8, ncol = 2)
dev.off()



