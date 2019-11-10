library(SNPRelate)
library(extrafont)
library(dendextend)
library(circlize)
library(pvclust)

## loading the gds of the data and pullling some attributes out
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data.gds")
# pop.code <- read.gdsn(index.gdsn(wheat.data, "MC"))
sample.id <- read.gdsn(index.gdsn(wheat.data, "sample.id"))

## pop.code from dendrogram clustering
load("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
dend <- genotypes.pvclust %>% as.dendrogram
pop.code <- vector()
pop.code[labels(dend)] <- c(rep("Red Spring", 238), rep("Chinese", 10), rep("Australian", 4),
                            rep("Extra Strong/\nPrairie", 28), rep("White Spring", 47), 
                            rep("Red Winter -\nOut", 3), rep("Red Winter", 58))
pop.code <- factor(pop.code, levels = c("Red Spring", "Chinese", "Australian", 
                                        "Extra Strong/\nPrairie", "White Spring", 
                                        "Red Winter -\nOut", "Red Winter"))

## LD-based SNP pruning
set.seed(1000)
snpset <- snpgdsLDpruning(wheat.data, ld.threshold=0.3, autosome.only = F,
                          maf = 0.05, missing.rate = 0.1, slide.max.bp = 345)
snpset.id <- unlist(snpset)

## IBS image
ibs <- snpgdsIBS(wheat.data, snp.id=snpset.id)
image(ibs$ibs[labels(dend), labels(dend)], col=terrain.colors(16))

## colour pallete for MC pop groups
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
palette(colours)

## MDS
loc <- cmdscale(1 - ibs$ibs, k = 4)
pairs(loc, col=pop.code, pch = 19, main = "Multidimensional Scaling Analysis (IBS)")

##hclustering on IBS
set.seed(1000)
ibs.hc <- snpgdsHCluster(ibs)
plot(ibs.hc$dendrogram)

hclust <- hclust(as.dist(1 - ibs$ibs), method = "average")
plot(hclust)
rv <- snpgdsCutTree(ibs.hc, z.threshold = 34)

## printing out the dendrogram
pdf("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Figures\\dend\\dend_test.pdf", family="Impact")
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 1)
circos.initialize("foo", xlim = c(0, 388), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:388, rep(0, 388), as.character(sample.id)[labels(dend)],
              facing = "clockwise", niceFacing = T, cex = 0.3, adj = c(0, -0.1), font = 2,
              col = sapply(as.numeric(pop.code), function(x) { return(colours[x]) } ))
}, track.height = 0.15, bg.border = NA)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:388, rep(0, 388), rv$sample.id[rv$samp.order],
              facing = "clockwise", niceFacing = T, cex = 0.3, adj = c(0, -0.1), font = 2,
              col = sapply(as.numeric(pop.code)[match(rv$samp.order, labels(dend))],
                           function(x) { return(colours[x]) } ))
}, track.height = 0.15, bg.border = NA)
max_height = attr(dend, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(rv$dendrogram, max_height = max_height)
}, track.height = 0.45, bg.border = NA)
legend("center", legend=c("Red Spring", "Chinese", "Australian", 
                          "Extra Strong/\nPrairie", "White Spring", 
                          "Red Winter -\nOut", "Red Winter"),
       pch=19, col=1:nlevels(pop.code), cex = 0.7)
text(0, 1, bquote(underline("Dendrogram of Wheats Coloured by Market class")))
circos.clear()
dev.off()

## IBD (kinda sucky)
ibd <- snpgdsIBDMoM(wheat.data, snp.id = snpset.id)
ibd.coeff <- snpgdsIBDSelection(ibd)
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

## Genetic Relationship Matrix
grm <- snpgdsGRM(wheat.data, snp.id = snpset.id, method = "Visscher", with.id = F)
rownames(grm) <- sample.id
snpgdsClose(wheat.data)
eig <- eigen(grm)
pairs(eig$vectors[,1:4], col=as.integer(pop.code), pch = 19, cex = 0.7, 
      main = "First Four PCs Plotted Against Each
      Other with Varieties Coloured by Market Class")



## probably not supposed to cluster on GRM
## hclust grm object
grm.hclust <- hclust(as.dist(3 - grm), method = "average")
plot(as.phylo(grm.hclust), type = "unrooted", show.tip.label = F)

## cicular plot
dend <- as.dendrogram(grm.hclust)
# pdf("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Figures\\dend\\dend_grm.pdf", family="Impact")
colours <- colors()[c(450, 139, 419, 53, 144, 554, 115, 26, 153, 49, 114)]
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 1)
circos.initialize("foo", xlim = c(0, 388), sector.width = 1)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:388, rep(0, 388), labels(dend),
              facing = "clockwise", niceFacing = T, cex = 0.3, adj = c(0, -0.1), font = 2,
              col = sapply(as.numeric(pop.code)[match(labels(dend), sample.id)], 
                           function(x) { return(colours[x]) } ))
}, track.height = 0.15, bg.border = NA)
dend_c <- color_branches(dend, k = 7, col = colors()[c(554, 144, 258, 53, 114, 450, 115)])
max_height = attr(dend_c, "height")
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(dend_c, max_height = max_height)
}, track.height = 0.45, bg.border = NA)
# legend("center", legend=levels(pop.code), pch=19, col=colours, cex = 0.45)
# text(0, 1, bquote(underline("Dendrogram of Wheats Coloured by Market class")))
circos.clear()
# dev.off()

## neighbourjoin grm object
library(ape)
install.packages("phangorn")
library(phangorn)
grm.nj <- nj(1 - grm)
pdf("C:\\Users\\Max_H\\OneDrive\\Pedagogy\\PGDA\\Figures\\dend\\dend_grm_nj.pdf", family="Impact")
plot(grm.nj, type = "fan", cex = 0.25)
dev.off()

