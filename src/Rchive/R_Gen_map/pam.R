library(SNPRelate)
# install.packages("fpc")
library(fpc)
library(extrafont)
library(circlize)
library(dendextend)
library(Hmisc)

wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data.gds")
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
genotypes <- read.gdsn(index.gdsn(wheat.data, "genotype"))
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\snp_indices.RData") 

## kmeans clustering
set.seed(1000)
# wss <- (nrow(t(genotypes[kept.indices,]))-1)*sum(apply(t(genotypes[kept.indices,]),2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(t(genotypes[kept.indices,]), 
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
     # ylab="Within groups sum of squares")
kmeans <- kmeans(t(genotypes[kept.indices,]), centers=2)

pamk <- pamk(dist(1-abs(cor(genotypes[kept.indices,]))), krange=2:15, diss = T)
# pamk <- pamk(t(genotypes[kept.indices,]), krange=2:6)
# pamk$pamobject$clustering

## pop.code from dendrogram clustering
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
dend <- genotypes.pvclust %>% as.dendrogram
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]

## loading the gds of the data and pullling some attributes out
sample.id <- read.gdsn(index.gdsn(wheat.data, "sample.id"))

## printing out the dendrogram
pdf("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Figures\\dend\\dend_pam.pdf", family="Impact")
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, 388), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:388, rep(0, 388), as.character(sample.id)[labels(dend)],
              facing = "clockwise", niceFacing = T, cex = 0.25, adj = c(0, -0.1), font = 2,
              col = sapply(pamk$pamobject$clustering[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.15, bg.border = NA)
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
palette(colours)
dend_colour <- color_branches(dend, k = 7, col = colours)
max_height <- attr(dend_colour, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(dend_colour, max_height = max_height)
}, track.height = 0.40, bg.border = NA)
legend("center", legend=c("Red Spring", "Chinese", "Australian", 
                          "Extra Strong/\nPrairie", "White Spring", 
                          "Red Winter -\nOut", "Red Winter"), 
       pch=19, col=1:7, cex = 0.4)
text(0, 1, bquote(underline("Dendrogram of Wheats Cut at 5 Highest Nodes")))
circos.clear()
dev.off()

## printing out the dendrogram
pdf("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Figures\\dend\\dend_kmeans.pdf", family="Impact")
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005))
circos.initialize("foo", xlim = c(0, 388), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:388, rep(0, 388), as.character(sample.id)[labels(dend)],
              facing = "clockwise", niceFacing = T, cex = 0.25, adj = c(0, -0.1), font = 2,
              col = sapply(kmeans$cluster[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.15, bg.border = NA)
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
palette(colours)
dend_colour <- color_branches(dend, k = 7, col = colours)
max_height <- attr(dend_colour, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(dend_colour, max_height = max_height)
}, track.height = 0.40, bg.border = NA)
legend("center", legend=c("Red Spring", "Chinese", "Australian", 
                          "Extra Strong/\nPrairie", "White Spring", 
                          "Red Winter -\nOut", "Red Winter"), 
       pch=19, col=1:7, cex = 0.4)
text(0, 1, bquote(underline("Dendrogram of Wheats Cut at 5 Highest Nodes")))
circos.clear()
dev.off()

## LD-based SNP pruning
set.seed(1000)
snpset <- snpgdsLDpruning(wheat.data, ld.threshold=0.3, autosome.only = F,
                          maf = 0.05, missing.rate = 0.1, slide.max.bp = 345)
snpset.id <- unlist(snpset)

## PCA of the different SNP sets and the full set
pca <- snpgdsPCA(wheat.data, sample.id = sample.id, snp.id = snpset.id,
                 num.thread = 4, autosome.only = F, maf = 0.05, missing.rate = 0.1)
pc.percent <- pca$varprop*100
snpgdsClose(wheat.data)

## colour pallete for oam pop groups
colours <- colors()[c(554, 144)]
palette(colours)

## Plotting the first four dimenions against each other
pdf("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Figures\\pca\\pca_pam.pdf", family="Impact")
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=pamk$pamobject$clustering, labels=lbls, pch = 19, cex = 0.45, 
      main = "First Four PCs Plotted Against Each Other
      With Varieties Coloured by 2 PAM Groups", xpd = F)
dev.off()

## colour pallete for kmeans pop groups
colours <- colors()[c(554, 144)]
palette(colours)

pdf("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Figures\\pca\\pca_kmeans.pdf", family="Impact")
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=kmeans$cluster, labels=lbls, pch = 19, cex = 0.45, 
      main = "First Four PCs Plotted Against Each Other
      With Varieties Coloured by 5 Dend Groups", oma = c(8,2,7,13), xpd = F)
dev.off()

snpgdsClose(wheat.data)