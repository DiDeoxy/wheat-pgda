library(extrafont)
library(SNPRelate)
library(plyr)

## useful data
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_subset_sample.gds")
snp.chrom <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snpgdsClose(wheat)
genes <- read.csv("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\final_genes_selected.csv",
                  header = F, stringsAsFactors = F, col.names = c("Name", "Contig", "Pos"))[,c(2,3,1)]
genes <- genes[order(genes[,1]),]
row.names(genes) <- NULL

replace <- 1:42
chromes <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
contigs <- as.vector(t(outer(chromes, c("part1", "part2"), paste, sep="_")))
names(replace) <- contigs

genes$Contig <- revalue(genes$Contig, replace)

genes_contigs <- data.frame()
count <- 1
number <- 0
for (i in 1:length(snp.chrom)) {
  if (count == snp.chrom[i]) {
    hits <- genes[which(genes$Contig == count),]
    number <- nrow(hits)
    if (nrow(hits) == 0) {
      temp <- data.frame(NA, NA, NA)
      names(temp) <- names(genes)
      genes_contigs <- rbind(genes_contigs, temp)
    } else {
      for (j in 0:dim(hits)[1]) {
        genes_contigs <- rbind(genes_contigs, hits[j,])
      }
    }
    if (number > 0) { number <- number - 1 }
    count = count + 1
  } else {
    if (number == 0) {
      temp <- data.frame(NA, NA, NA)
      names(temp) <- names(genes)
      genes_contigs <- rbind(genes_contigs, temp)
    } else {
      number <- number - 1
    }
  }
}

## Plotting phi by positons and chr
# some values that are reused on the graphing function
cex <- 0.8
pch <- 20
ylim <- c(0, 1)
num.snps <- table(snp.chrom)
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
order <- c(rep(1,6), rep(2,6), rep(3,6), rep(4,6), rep(5,6), rep(6,6), rep(7,6))
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
plot.wheat <- function (chr_markers, i) { plot(chr_markers, pch = pch, ylim = ylim, 
                                               xlim = c(0, max(snp.pos)), cex.axis = cex, 
                                               col = colours[order[i]], yaxt= "n", xaxt = "n") }

## loading pre-calcualted per locus amova results
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\amova_hrs_hrw.RData")

## finding phi resulsts for each locus and the top quantiles
phi.all <- vector()
for (i in 1:length(wheat.amova)) {
  if (is.null(wheat.amova[[i]]$statphi[1,1])) {
    phi.all[i] <- NA
  } else {
    phi.all[i] <- abs(wheat.amova[[i]]$statphi[1,1])
  }
}
signif <- quantile(phi.all, prob = 0.95, na.rm = T)
signif2 <- quantile(phi.all, prob = 0.99, na.rm = T)

## plotting
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\loci\\HRS_amova_hrs_vs_hrw.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
par(mfrow = c(7,6), oma = c(3,3,4,2))
# counts how many times tapply is called
count <- 1
by(cbind(snp.pos, phi.all, genes_contigs), snp.chrom, function (chr_markers) {
  if (count %in% seq(1, 5, 2)) { # the first graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot.wheat(chr_markers[,1], chr_markers[,2], count)
    abline(chr_markers[,4], col = "black")
    mtext(text = "Long", 3, line = 0.5, cex = cex)
    mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    if (count == 5) {
      text(genes[1,2], labels = genes[1,3])
      
    }
  } else if (count %in% seq(2, 6, 2)) { # the first graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    plot.wheat(chr_markers, count)
    axis(4, at=seq(0, 1, 0.1))
    mtext(text = "Short", 3, line = 0.5, cex = cex)
  } else if (count %in% seq(7, 41, 2)) { # the remaining graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot.wheat(chr_markers, count)
    mtext(text = labels[count - (count - 1)/2], 2, line = 0.5, cex = cex)
    if (count == 7) {
      text(genes[2,2], labels = genes[2,3])
      abline(v = genes[2,2], col = "black")
    }
    if (count == 9) {
      text(genes[3,2], labels = genes[3,3])
      abline(v = genes[3,2], col = "black")
    }
    if (count == 11) {
      text(genes[4,2], labels = genes[4,3])
      abline(v = genes[4,2], col = "black")
    }
    if (count == 13) {
      text(genes[9,2], labels = genes[9,3])
      abline(v = genes[9,2], col = "black")
    }
    if (count == 15) {
      # text(genes[9,2], labels = genes[9,3])
      abline(v = genes[11,2], col = "black")
      abline(v = genes[12,2], col = "black")
      abline(v = genes[13,2], col = "black")
    }
    if (count == 17) {
      text(genes[15,2], labels = genes[9,3])
      abline(v = genes[15,2], col = "black")
    }
    if (count == 21) {
      abline(v = genes[18,2], col = "black")
    }
    if (count == 23) {
      abline(v = genes[19,2], col = "black")
    }
    if (count == 29) {
      abline(v = genes[22,2], col = "black")
      abline(v = genes[23,2], col = "black")
    }
    if (count == 39) {
      abline(v = genes[25,2], col = "black")
    }
  } else { # the remaining graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    plot.wheat(chr_markers, count)
    if (count %in% c(14, 16, 18, 26, 28, 30, 38, 40, 42)) axis(4, at=seq(0, 1, 0.1))
    if (count == 12) {
      # text(genes[2,2], labels = genes[2,3])
      abline(v = genes[5,2], col = "black")
      abline(v = genes[6,2], col = "black")
      abline(v = genes[7,2], col = "black")
      abline(v = genes[8,2], col = "black")
    }
    if (count == 14) {
      text(genes[10,2], labels = genes[10,3])
      abline(v = genes[10,2], col = "black")
    }
    if (count == 16) {
      text(genes[14,2], labels = genes[14,3])
      abline(v = genes[14,2], col = "black")
    }
    if (count == 18) {
      text(genes[16,2], labels = genes[16,3])
      abline(v = genes[16,2], col = "black")
    }
    if (count == 20) {
      abline(v = genes[17,2], col = "black")
    }
    if (count == 26) {
      abline(v = genes[20,2], col = "black")
    }
    if (count == 28) {
      abline(v = genes[21,2], col = "black")
    }
    if (count == 30) {
      abline(v = genes[24,2], col = "black")
    }
  }
  abline(h = signif, col = "blue")
  abline(h = signif2, col = "red")
  # mtext(3, text = paste(labels[count], ": ", num.snps[count]),
  #       line = -1.25, cex = cex)
  count <<- count + 1
  # which(x[,2] >= signif2)
})
title(xlab = "Marker by Chromosome Position", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5)
title("AMOVA Phi Statistic By Locus Between HRS and HRW Wheats", outer = T, cex.main = 1.5)
dev.off()

