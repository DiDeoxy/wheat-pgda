library(SNPRelate)
library(dendextend)
library(ape)
library(plyr)
# library(fpc)
library(extrafont)
# library(dendextend)

## loading the gds of the data and pullling some attributes out, etc.
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_subset_sample.gds")
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
genotypes <- replace(genotypes, genotypes == 3, NA)


## Getting the nj order
neighbour <- nj(1 - snpgdsIBS(wheat, autosome.only = F)$ibs)
snpgdsClose(wheat)
labels(neighbour) <- sample.id
neighbour <- root(neighbour, which(sample.id == "Selkirk"), r = TRUE)
neighbour <- chronos(neighbour, lambda = 0, model = "correlated")
neighbour <- as.dendrogram(neighbour)

## getting the nj pop codes
pop.code <- vector()
pop.code[match(labels(neighbour), sample.id)] <- c(rep("Hard Red Spring 5", 83), rep("Hard Red Spring 1-4", 119),
                                                   rep("Prairie Spring Red", 32), rep("Out", 1))
pop.code <- factor(pop.code)
sample.compare <- c(which(pop.code == "Hard Red Spring 5"), which(pop.code == "Hard Red Spring 1-4"))

wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_subset_sample.gds")
## calculating fst in sliding window along whole genome
fst <- snpgdsSlidingWindow(wheat, FUN = "snpgdsFst", unit = "locus", winsize = 1, shift = 1, 
                           population = factor(pop.code[sample.compare]), 
                           sample.id = sample.id[sample.compare],
                           method ="W&H02", win.start = 0, as.is = "numeric")
snpgdsClose(wheat)

## finding fst quantiles
fst_all <- vector()
for (i in c(seq(3, 167, 4))) {
  fst_all <- c(fst_all, fst[i])
}
fst_all <- unlist(fst_all)
signif <- quantile(fst_all, prob = 0.95)
signif2 <- quantile(fst_all, prob = 0.99)

## Some constant values for the plot
ylim <- c(0,1)
cex <- 0.8
pch <- 20
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
order <- c(rep(1,6), rep(2,6), rep(3,6), rep(4,6), rep(5,6), rep(6,6), rep(7,6))
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]

## plotting
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\loci\\HRS_fst_HRS1-4_vs_HRS5.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
par(mfrow = c(7,6), oma = c(3,3,4,2))
count <- 1
for (i in seq(3, 167, 4)) {
  if (i %in% c(3, 11, 19)) { # the first graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot(unlist(fst[i+2]), unlist(fst[i]), ylim = ylim, cex.axis = cex, pch = pch,
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    mtext(text = "Long", 3, line = 0.5, cex = cex)
    mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
  } else if (i %in% c(7, 15, 23)) { # the first graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    plot(unlist(fst[i+2]), unlist(fst[i]), ylim = ylim, cex.axis = cex, pch = pch, 
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    axis(4, at=seq(0, 1, 0.1))
    mtext(text = "Short", 3, line = 0.5, cex = cex)
  } else if (i %in% seq(27, 163, 8)) { # the remaining graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot(unlist(fst[i+2]), unlist(fst[i]), ylim = ylim, cex.axis = cex, pch = pch,
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
  } else {
    par(mar = c(0, 0, 0, 2))
    plot(unlist(fst[i+2]), unlist(fst[i]), ylim = ylim, cex.axis = cex, pch = pch, 
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    if (count %in% c(14, 16, 18, 26, 28, 30, 38, 40, 42)) axis(4, at=seq(0, 1, 0.1))
  }
  abline(h = signif, col = "blue")
  abline(h = signif2, col = "red")
  count <- count + 1
}
title(xlab = "Marker by Position in BPs", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "Wright's Fst Value", outer = T, cex.lab = 1.5, line=0.5)
title("Fst by Chromosome Between HRS1-4 and HRS5", outer = T, cex.main = 1.5)
dev.off()

