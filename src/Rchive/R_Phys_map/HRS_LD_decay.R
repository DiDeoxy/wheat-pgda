# install.packages("sommer")
library(sommer)
library(SNPRelate)
library(extrafont)

## loading the gds of the data and pullling some attributes out
wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\HRS_phys_subset_sample.gds")
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
snp.id <- read.gdsn(index.gdsn(wheat, "snp.id"))
snp.chrom <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snpgdsClose(wheat)

genotypes <- replace(genotypes, genotypes == 0, -1)
genotypes <- replace(genotypes, genotypes == 2, 1)
genotypes <- replace(genotypes, genotypes == 3, NA)
rownames(genotypes) <- snp.id
colnames(genotypes) <- sample.id

map <- data.frame(Locus = snp.id, LG = snp.chrom, Position = snp.pos)

# wheat.decay <- LD.decay(t(genotypes), map)
# wheat.decay.unlinked <- LD.decay(t(genotypes), map, unlinked = TRUE, gamma = .95)


half.decay.distance <- list()
plot.wheat <- function () {
    with(wheat.decay$by.LG[[i]][which(wheat.decay$by.LG[[i]]$p < .001),],
         plot(r2~d,col=colours[order[i]],
              xlim=c(0,max(wheat.decay$by.LG[[i]]$d)), ylim=c(0,1),
              pch=20, cex=0.5, yaxt="n",
              xaxt = "n", ylab=expression(r^2),
              xlab="Distance in BPs"))
    values <- wheat.decay$by.LG[[i]][which(wheat.decay$by.LG[[i]]$p < .001),]
    mod <- loess(r2 ~ d, data = values)
    par(new=T)
    xl <- seq(min(values$d), max(values$d), (max(values$d) - min(values$d))/1000)
    lilo <- predict(mod)
    lines(xl, predict(mod, xl), col="green", lwd=2, bty="n", xaxt="n", yaxt="n",
          type="l", ylim=c(0,1), ylab="", xlab="")
    abline(h=wheat.decay.unlinked$by.LG[[i]], col="blue")
}

png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\loci\\HRS_LD_decay.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
cex <- 0.8
pch <- 20
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
order <- c(rep(1,6), rep(2,6), rep(3,6), rep(4,6), rep(5,6), rep(6,6), rep(7,6))
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
par(mfrow = c(7,6), oma = c(3,3,4,2))
for (i in 1:42) {
  if (i %in% seq(1, 5, 2)) { # the first graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot.wheat()
    mtext(text = "Long", 3, line = 0.5, cex = cex)
    mtext(text = labels[i - (i - 1)/2], 2, line = 0.5, cex = cex)
  } else if (i %in% seq(2, 6, 2)) { # the first graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    if (i == 6){
      plot(0,xaxt='n',yaxt='n',pch='',ylab='',xlab='', ylim = c(0,1))
    } else {
      plot.wheat()
    }
    axis(4, at=seq(0, 1, 0.1))
    mtext(text = "Short", 3, line = 0.5, cex = cex)
  } else if (i %in% seq(7, 41, 2)) { # the remaining graphs of coulmns 1, 3, and 5
    par(mar = c(0, 2, 0, 0))
    plot.wheat()
    mtext(text = labels[i - (i - 1)/2], 2, line = 0.5, cex = cex)
  } else { # the remaining graphs of coulmns 2, 4, and 6
    par(mar = c(0, 0, 0, 2))
    plot.wheat()
    if (i %in% c(14, 16, 18, 26, 28, 30, 38, 40, 42)) axis(4, at=seq(0, 1, 0.1))
  }
}
title(xlab = "Physical Distance Between Markers", outer = T, cex.lab = 1.5, line = 1)
title(ylab = "Level of LD", outer = T, cex.lab = 1.5, line = 0.5)
title("Distribution of LD Across Chromosomes", outer = T, cex.main = 2.5, line = 2)
dev.off()

1# mean(unlist(half.decay.distance), na.rm = T)

r2.ave <- list()
for (i in 1:42) {
  if (i %in% c(1, 5, 6, 7, 12, 15, 17, 18, 19, 23, 24, 25, 29, 30, 31, 32, 35, 36, 42)) {
    next
  } else {
    r2.ave[[i]] <- mean(wheat.decay$by.LG[[i]][which(wheat.decay$by.LG[[i]]$p < .001),]$r2)
  }
}
r2.ave
min(unlist(r2.ave))
max(unlist(r2.ave))
mean(unlist(r2.ave))
