library(SNPRelate)
library(fpc)
library(extrafont)
library(devEMF)

## loading the gds of the data and pullling some attributes out, etc.
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data_2.gds")
genotypes <- read.gdsn(index.gdsn(wheat.data, "genotype"))
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
pop.code <- read.gdsn(index.gdsn(wheat.data, "season"))
sample.id <- read.gdsn(index.gdsn(wheat.data, "sample.id"))

# Removing SNPs with high MR and low MAF to create informative snp set
snpset <- snpgdsSelectSNP(wheat.data, autosome.only = F, maf = 0.05, missing.rate = 0.1)
snp.set.id <- unlist(snpset)
informative <- match(snp.set.id, snp.id)

## reordered chromosomes for snp positions of informative snp set (for plotting genomes together)
snp.pos <- as.numeric(read.gdsn(index.gdsn(wheat.data, "snp.position")))[informative] / 1000
snp.chr <- factor(read.gdsn(index.gdsn(wheat.data, "snp.chromosome"))[informative], 
                  levels = c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3)))
snp.pos.2 <- tapply(snp.pos, snp.chr, function (x) {
  x
}, simplify = F)

## creating partition groups using pruned snp set
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\snp_indices.RData")
pamk <- pamk(dist(1-abs(cor(genotypes[kept.indices,]))), krange=2:15, diss = T)

## calculating fst in sliding window along whole genome using informative snp set
fst <- snpgdsSlidingWindow(wheat.data, FUN = "snpgdsFst", unit = "locus", winsize = 1, shift = 1, 
                           sample.id = as.character(sample.id[-c(which(pop.code == "UNKNOWN"), which(pop.code == "N/A"))]),
                           population = factor(pop.code[-c(which(pop.code == "UNKNOWN"), which(pop.code == "N/A"))]),
                           snp.id = snp.set.id, method ="W&H02", win.start = 0, as.is = "numeric")
snpgdsClose(wheat.data)

## finding fst quantiles
fst_all <- vector()
for (i in c(seq(3, 75, 12),seq(7, 79, 12),seq(11, 83, 12))) {
  fst_all <- c(fst_all, fst[[i]])
}
signif <- quantile(fst_all, prob = 0.95)
signif2 <- quantile(fst_all, prob = 0.99)

## Some constant values for the plot
ylim <- c(0,1)
cex <- 0.8
pch <- 20
labels <- c("NA", "1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", 
            "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", 
            "7B", "7D")[c(1, seq(2, 20, 3), seq(3, 21, 3), seq(4, 22, 3))]
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]

## plotting
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\loci\\fst\\fst_full_season.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
par(mfrow = c(3,7), oma = c(3,5,3,2), mar = c(3, 0, 0, 0))
count = 1
pos = 0
for (i in c(seq(3,75,12),seq(7,79,12),seq(11,83,12))) {
  if (i %in% c(3,7,11)) {
    plot(snp.pos.2[[count]][1:length(fst[[i]])], fst[[i]], ylim = ylim, cex.axis = cex, pch = pch,
         col = colours[1 + count %% 7])
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
  } else {
    plot(snp.pos.2[[count]][1:length(fst[[i]])], fst[[i]], ylim = ylim, cex.axis = cex, pch = pch,
         col = colours[1 + count %% 7], yaxt = "n")
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
  }
  mtext(3, text = paste(labels[count + 1], ": ", length(fst[[i]])),
        line = -1.25, cex = cex)
  count <- count + 1
  pos <- pos + length(fst[[i]])
}
title(xlab = "Marker by Position in cM", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "Wright's Fst Value", outer = T, cex.lab = 1.5)
title("Figure 7: Fst by Marker Position Between Season Type", outer = T, cex.main = 1.5)
dev.off()