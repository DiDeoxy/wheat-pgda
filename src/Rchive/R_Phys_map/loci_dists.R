library(SNPRelate)
library(extrafont)

wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_subset_both.gds")
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chrom <- factor(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
snpgdsClose(wheat)

num.snps <- table(snp.chrom)
xlim = c(0, max(snp.pos))
cex <- 0.8
pch <- 20
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
order <- c(rep(1,6), rep(2,6), rep(3,6), rep(4,6), rep(5,6), rep(6,6), rep(7,6))
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\loci\\HRS_pos.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
par(mfrow = c(7,6), oma = c(5,3,4,2))
count <- 1
by(snp.pos, snp.chrom, function (chr) {
  if (count %in% seq(1, 5, 2)) { # the first graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot(x = chr, y = 1:length(chr), pch = "|", cex.axis = cex, xlim = xlim,
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    mtext(text = "Long", 3, line = 0.5, cex = cex)
    mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    mtext(text = paste("SNPs: ", num.snps[count]), 1, line = -1, cex = 0.7, adj = 1)
  } else if (count %in% seq(2, 6, 2)) { # the first graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    plot(x = chr, y = 1:length(chr), pch = "|", cex.axis = cex,  xlim = xlim,
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    mtext(text = "Short", 3, line = 0.5, cex = cex)
    mtext(text = paste("SNPs: ", num.snps[count]), 1, line = -1, cex = 0.7, adj = 1)
  } else if (count %in% seq(7, 41, 2)) { # the remaining graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot(x = chr, y = 1:length(chr), pch = "|", cex.axis = cex,  xlim = xlim,
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    mtext(text = paste("SNPs: ", num.snps[count]), 1, line = -1, cex = 0.7, adj = 1)
    if (count %in% c(37, 39, 41)) axis(1, las = 2, cex.axis = 0.75)
  } else { # the remaining graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    plot(x = chr, y = 1:length(chr), pch = "|", cex.axis = cex,  xlim = xlim,
         col = colours[order[count]], xaxt = "n", yaxt = "n")
    mtext(text = paste("SNPs: ", num.snps[count]), 1, line = -1, cex = 0.7, adj = 1)
    if (count %in% c(38, 40, 42)) axis(1, las = 2, cex.axis = 0.75)
  }
  count <<- count + 1
})
title(xlab = "Marker Position in Base Pairs", outer = T, cex.lab = 1.5, line = 3.5)
title(ylab = "Index Position", outer = T, cex.lab = 1.5, line = 0.5)
title("Disrtribution of Pruned Markers Across Chromosomes", outer = T, cex.main = 2.5, line = 2)
dev.off()
