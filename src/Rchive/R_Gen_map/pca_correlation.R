library(SNPRelate)
library(extrafont)
library(dendextend)
library(rgl)
library(plyr) #using revalue function
library(devEMF)


## loading the gds of the data and pullling some attributes out
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data_2.gds")

## pop.code from dendrogram clustering
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\genotypes_pruned_pvclust.RDATA")
dend <- genotypes.pvclust %>% as.dendrogram
pop.code <- vector()
pop.code[labels(dend)] <- c(rep("Red Spring", 238), rep("Chinese", 10), rep("Australian", 4),
                            rep("Extra Strong/\nPrairie", 28), rep("White Spring", 47), 
                            rep("Red Winter-\nOut", 3), rep("Red Winter", 58))
pop.code <- factor(pop.code, levels = c("Red Spring", "Chinese", "Australian", 
                                        "Extra Strong/\nPrairie", "White Spring", 
                                        "Red Winter-\nOut", "Red Winter"))
snp.chr <- factor(read.gdsn(index.gdsn(wheat.data, "snp.chromosome")),
                  levels = c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3)))
snp.pos <- read.gdsn(index.gdsn(wheat.data, "snp.position"))/1000
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
snpset <- snpgdsSelectSNP(wheat.data, autosome.only = F, maf = 0.05, missing.rate = 0.1)
snp.set.id <- unlist(snpset)
informative <- match(snp.set.id, snp.id)
snp.chr.inf <- snp.chr[informative]
snp.pos.inf <- snp.pos[informative]
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\snp_indices.RData")

## LD-based SNP pruning
set.seed(1000)
snpset <- snpgdsLDpruning(wheat.data, ld.threshold=0.3, autosome.only = F,
                          maf = 0.05, missing.rate = 0.1, slide.max.bp = 345)
snpset.id <- unlist(snpset)

## PCA corr of the SNP set
pca <- snpgdsPCA(wheat.data, snp.id = snpset.id, num.thread = 4, autosome.only = F, 
                 maf = 0.05, missing.rate = 0.1)
pc.percent <- pca$varprop*100
corr <- snpgdsPCACorr(pca, wheat.data, eig.which=1:3, snp.id = snp.set.id)
snpgdsClose(wheat.data)

## look at the correlation of first 3 PCs with SNPs
# some common variables
num.snps <- table(snp.chr.inf)
cex <- 0.8
ylim = c(0,1)
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
# Chromosome labels
labels <- c("NA", "1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", 
            "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D")
# genome chromosome label order
labels.reordered <- labels[c(seq(2, 20, 3), seq(3, 21, 3), seq(4, 22, 3), 1)]

# making the graph itself
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\loci\\pca\\correlation_pc1.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
## PC1
pc <- 1
## significant quantiles
signif <- quantile(corr$snpcorr[pc,], prob = 0.95)
signif2 <- quantile(corr$snpcorr[pc,], prob = 0.99)
# making 21 graphs
par(mfrow = c(3,7), oma = c(3,5,3,2), mar = c(3, 0, 0, 0))
# counts how many times tapply is called
count = 1
## plotting
foo <- by(cbind(snp.pos.inf, abs(corr$snpcorr[pc,])), snp.chr.inf, function (x) {
  # the first graph in each row
  if (count %% 7 == 1) {
    plot(x, pch = 20, las = 3, 
         ylim = ylim,  col = colours[1 + count %% 7])
    mtext(3, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = -1.25, cex = cex)
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
  }
  # the last 6 graphs in each row
  if (count %% 7 != 1) {
    plot(x, pch = 20, las = 3, yaxt = "n", 
         ylim = ylim, xlab = labels.reordered[count],  col = colours[1 + count %% 7])
    mtext(3, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = -1.25, cex = cex)
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
  }
  count <<- count + 1
})
# titles and axes labels
title("Figure 5: SNP Correlation By Marker Position to PC1", out = T, cex.main = 1.5)
title(xlab = "Marker by Position in cM", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "Absolute Correlation Coefficient", outer = T, cex.lab = 1.5)
dev.off()

## PC2
pc <- 2
signif <- quantile(corr$snpcorr[pc,], prob = 0.95)
signif2 <- quantile(corr$snpcorr[pc,], prob = 0.99)
# making the graph itself
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\loci\\pca\\correlation_pc2.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
# making 21 graphs
par(mfrow = c(3,7), oma = c(3,5,2,1))
# counts how many times tapply is called
count = 1
## plotting
foo <- by(cbind(snp.pos.inf, abs(corr$snpcorr[pc,])), snp.chr.inf, function (x) {
  par(mar = c(2, 0, 0, 0))
  # the first graph in each row
  if (count %% 7 == 1) {
    plot(x, pch = 20, las = 3, xaxt = "n", 
         ylim = ylim,  col = colours[1 + count %% 7])
    mtext(1, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = 0.5, cex = cex)
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
  }
  # the last 6 graphs in each row
  if (count %% 7 != 1) {
    plot(x, pch = 20, las = 3, xaxt = "n", yaxt = "n", 
         ylim = ylim, xlab = labels.reordered[count],  col = colours[1 + count %% 7])
    mtext(1, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = 0.5, cex = cex)
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
  }
  count <<- count + 1
})
# titles and axes labels
title("PCA Correlation By Locus to PC2", out = T, cex.main = 1.5)
title(xlab = "Marker by Position in cM", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "Absolute Correlation Coefficient", outer = T, cex.lab = 1.5)
dev.off()

## PC3
pc <- 3
signif <- quantile(corr$snpcorr[pc,], prob = 0.95)
signif2 <- quantile(corr$snpcorr[pc,], prob = 0.99)
# making the graph itself
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\loci\\pca\\correlation_pc3.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
# making 21 graphs
par(mfrow = c(3,7), oma = c(3,5,2,1))
# counts how many times tapply is called
count = 1
## plotting
foo <- by(cbind(snp.pos.inf, abs(corr$snpcorr[pc,])), snp.chr.inf, function (x) {
  par(mar = c(2, 0, 0, 0))
  # the first graph in each row
  if (count %% 7 == 1) {
    plot(x, pch = 20, las = 3, xaxt = "n", 
         ylim = ylim,  col = colours[1 + count %% 7])
    mtext(1, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = 0.5, cex = cex)
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
  }
  # the last 6 graphs in each row
  if (count %% 7 != 1) {
    plot(x, pch = 20, las = 3, xaxt = "n", yaxt = "n", 
         ylim = ylim, xlab = labels.reordered[count],  col = colours[1 + count %% 7])
    mtext(1, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = 0.5, cex = cex)
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
  }
  count <<- count + 1
})
# titles and axes labels
title("PCA Correlation By Locus to PC3", out = T, cex.main = 1.5)
title(xlab = "Marker by Position in cM", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "Absolute Correlation Coefficient", outer = T, cex.lab = 1.5)
dev.off()