library(SNPRelate)
library(extrafont)
library(dendextend)
# library(rgl)
library(plyr) #using revalue function
library(ape)
# install.packages("GeneticSubsetter")
library(GeneticSubsetter)
library(lintr)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_imputed_subset_sample.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chrom <- factor(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
snpgdsClose(wheat)

genotypes <- replace(genotypes, genotypes == 0, -1)
genotypes <- replace(genotypes, genotypes == 2, 1)
genotypes <- replace(genotypes, genotypes == 3, 0)

PICs <- vector()
for (i in 1:dim(genotypes)[1]) {
  PICs <- c(PICs, PicCalc(t(genotypes[i,])))
}


ylim <- c(0,0.5)
cex <- 0.8
pch <- 20
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
order <- c(rep(1,6), rep(2,6), rep(3,6), rep(4,6), rep(5,6), rep(6,6), rep(7,6))
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\loci\\PIC\\HRS_PIC.png", 
    family="Calibri", width = 6, height = 6, pointsize = 5, units = "in", res = 300)
par(mfrow = c(7,6), oma = c(3,5,3,2))
count <- 1
by(cbind(snp.pos, PICs), snp.chrom, function (chr) {
  if (count %in% seq(1, 5, 2)) { # the first graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot(chr, pch = 20, ylim = ylim, cex.axis = cex, 
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    mtext(text = "Part 1", 3, line = 0.5, cex = cex)
    if (count %% 2 == 1) {
      mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    }
  } else if (count %in% seq(2, 6, 2)) { # the first graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    plot(chr, pch = 20, ylim = ylim, cex.axis = cex, 
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    axis(4, at=seq(0, 0.5, 0.1))
    mtext(text = "Part 2", 3, line = 0.5, cex = cex)
  } else if (count %in% seq(7, 41, 2)) { # the remaining graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot(chr, pch = 20, ylim = ylim, cex.axis = cex, 
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    if (count %% 2 == 1) {
      mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    }
  } else { # the remaining graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    plot(chr, pch = 20, ylim = ylim, cex.axis = cex, 
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    axis(4, at=seq(0, 0.5, 0.1))
    if (count %% 2 == 1) {
      mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    }
  }
  count <<- count + 1
})
dev.off()
