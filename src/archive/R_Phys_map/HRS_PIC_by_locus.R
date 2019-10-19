library(SNPRelate)
library(extrafont)
library(dendextend)
# library(rgl)
library(plyr) #using revalue function
library(ape)
# install.packages("GeneticSubsetter")
library(GeneticSubsetter)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_imputed_subset_sample.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snp.chrom <- factor(read.gdsn(index.gdsn(wheat, "snp.chromosome")),
                    c(seq(1, 37, 6), seq(2, 38, 6),
                      seq(3, 39, 6), seq(4, 40, 6),
                      seq(5, 41, 6), seq(6, 42, 6)
                      ))
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
colours <- rep(colors()[c(554, 144, 258, 53, 114, 450, 115)],6)
labels <- as.vector(t(outer(c("A", "B", "D"), c("part1", "part2"), paste, sep="_")))
labels <- as.vector(outer(as.character(1:7), labels, paste, sep=""))
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\loci\\PIC\\pic_hrs.png", 
    family="Calibri", width = 6, height = 6, pointsize = 5, units = "in", res = 300)
par(mfrow = c(6,7), oma = c(3,5,3,2))
count <- 1
count2 <- 1
by(PICs, snp.chrom, function (chr) {
  if (count %in% c(1,15,29)) { # the first graph of rows 1, 3, and 5
      par(mar = c(0, 0, 2, 0))
      plot(chr, pch = 20, ylim = ylim, cex.axis = cex, 
           col = colours[count], xaxt = "n")
      } else if (count %in% c(8,22,36)) { # the first graph of rows 2, 4, and 6
      par(mar = c(2, 0, 0, 0))
      plot(chr, pch = 20, ylim = ylim, cex.axis = cex, 
           col = colours[count], xaxt = "n")
      mtext(1, text = paste(labels[count2], ": ", length(chr)),
            line = 1, cex = cex)
      count2 <<- count2 + 1
      } else if (count %in% c(2:7,16:21,30:35)) { # the six graphs following the first in rows 1, 3, and 5
      par(mar = c(0, 0, 2, 0))
      plot(chr, pch = 20, ylim = ylim, cex.axis = cex, 
           col = colours[count], xaxt = "n", yaxt = "n")
      } else { # the six graphs following the first in rows 2, 4, and 6
      par(mar = c(2, 0, 0, 0))
      plot(chr, pch = 20, ylim = ylim, cex.axis = cex, 
           col = colours[count], xaxt = "n", yaxt = "n")
      mtext(1, text = paste(labels[count2], ": ", length(chr)),
            line = 1, cex = cex)
      count2 <<- count2 + 1
      }
  count <<- count + 1
})
dev.off()
