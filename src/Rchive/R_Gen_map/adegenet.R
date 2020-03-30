library(adegenet)
# install.packages("adegenet")
library(pvclust)
library(circlize)
library(extrafont)
library(dendextend)

## heirarchical clustering

## loading genlight object wheat.data
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genlight.RData")
pop.code <- wheat.data$other$mc

##loading hclust dendrogram
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
dend <- genotypes.pvclust %>% as.dendrogram

## creating the dendrogram
t_dend <- as.dendrogram(pvclust(t(as.matrix(wheat.data)), method.dist = "euclidean", method.hclust = "average", 
                 nboot = 1, r = 1))

## printing out the dendrogram
pdf("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Figures\\dend\\dend_nei.pdf", family="Impact")
# colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
colours <- colors()[c(450, 139, 419, 53, 144, 554, 115, 26, 153, 49, 114)]
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 1)
circos.initialize("foo", xlim = c(0, 388), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:388, rep(0, 388), labels(t_dend),
              facing = "clockwise", niceFacing = T, cex = 0.3, adj = c(0, 0), font = 2,
                col = sapply(as.numeric(pop.code)[match(labels(t_dend), indNames(wheat.data))], 
                             function(x) { return(colours[x]) } ))
  # col = sapply(as.numeric(pop.code)[labels(dend)], function(x) { return(colours[x]) } ))
}, track.height = 0.15, bg.border = NA)
max_height = attr(t_dend, "height")
t_dend = color_branches(t_dend, k = 7, col = colors()[c(554, 144, 258, 53, 114, 450, 115)])
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(t_dend, max_height = max_height)
}, track.height = 0.45, bg.border = NA)
legend("center", legend=levels(pop.code), pch=19, col=colours, cex = 0.45)
circos.clear()
dev.off()

## variations
# library(poppr)
# t_dend <- as.dendrogram(aboot(wheat.data, tree = "upgma", distance = "edwards.dist", showtree = F, sample = 1))

library(StAMPP)
# install.packages("StAMPP")
t_dist <- as.dist(stamppNeisD(wheat.data, pop = F))
t_dend <- as.dendrogram(hclust(t_dist, method = "average"))

##### other stuff
##
glPlot(wheat.data[labels(dend),], posi="topleft")

##
myFreq <- glMean(wheat.data)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)

myFreq <- glMean(wheat.data)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)

pca <- glPca(wheat.data, parallel = F)
scatter(pca)

## k-means cluster of PCA
grp <- find.clusters(wheat.data, max.n.clust=30, parallel = F, n.pca = 388)
table(pop.code, grp$grp)

## DPCA need to figure out the right number of clusters, gettin overlap like mad here
dapc1 <- dapc(wheat.data, grp$grp, parallel = F, n.pca = 6, n.da = 7)

op <- par(mfrow = c(3,3))
scatter(dapc1, 2, 1, scree.da=F, scree.pca = F, bg="white", posi.pca="bottomleft", legend=F,
        txt.leg=paste("group", 1:14),
        col=colors()[c(116, 554, 53, 26, 100, 12, 490, 115, 42, 139, 451, 142, 114, 195, 557, 79, 153, 48)])
scatter(dapc1, 3, 1, scree.da=F, scree.pca = F, bg="white", posi.pca="bottomleft", legend=F,
        txt.leg=paste("group", 1:14),
        col=colors()[c(116, 554, 53, 26, 100, 12, 490, 115, 42, 139, 451, 142, 114, 195, 557, 79, 153, 48)])
scatter(dapc1, 4, 1, scree.da=F, scree.pca = F, bg="white", posi.pca="bottomleft", legend=F,
        txt.leg=paste("group", 1:14),
        col=colors()[c(116, 554, 53, 26, 100, 12, 490, 115, 42, 139, 451, 142, 114, 195, 557, 79, 153, 48)])
plot.new()
scatter(dapc1, 3, 2, scree.da=F, scree.pca = F, bg="white", posi.pca="bottomleft", legend=F,
        txt.leg=paste("group", 1:14),
        col=colors()[c(116, 554, 53, 26, 100, 12, 490, 115, 42, 139, 451, 142, 114, 195, 557, 79, 153, 48)])
scatter(dapc1, 4, 2, scree.da=F, scree.pca = F, bg="white", posi.pca="bottomleft", legend=F,
        txt.leg=paste("group", 1:14),
        col=colors()[c(116, 554, 53, 26, 100, 12, 490, 115, 42, 139, 451, 142, 114, 195, 557, 79, 153, 48)])
plot.new()
plot.new()
scatter(dapc1, 4, 3, scree.da=F, scree.pca = F, bg="white", posi.pca="bottomleft", legend=F,
        txt.leg=paste("group", 1:14),
        col=colors()[c(116, 554, 53, 26, 100, 12, 490, 115, 42, 139, 451, 142, 114, 195, 557, 79, 153, 48)])
par(op)




