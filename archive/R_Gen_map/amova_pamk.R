library(poppr)
library(SNPRelate)
library(extrafont)
library(devEMF)

## loading pre-calcualted per locus amova results
load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\AMOVA\\amova_pk.RData")

## finding phi resulsts for each locus and the top quantiles
phi.all <- vector()
for (i in 1:length(wheat.amova)) {
  phi.all <- c(phi.all, wheat.amova[[i]]$statphi[1,1])
}
signif <- quantile(phi.all, prob = 0.95)
signif2 <- quantile(phi.all, prob = 0.99)

## useful data
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data_2.gds")
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
snp.pos <- snp.pos <- read.gdsn(index.gdsn(wheat.data, "snp.position"))/1000
snp.chrom <- read.gdsn(index.gdsn(wheat.data, "snp.chromosome"))
# finding informative subset of data
snpset <- snpgdsSelectSNP(wheat.data, autosome.only = F, maf = 0.05, missing.rate = 0.1)
snpset.id <- unlist(snpset)
informative <- match(snpset.id, snp.id)
snpgdsClose(wheat.data)
# subsetting data and reordering chromosomes for graphing
snp.chrom <- factor(snp.chrom[informative],levels = c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3)))
snp.pos <- snp.pos[informative]

## making data frame for graphing
pos.phi <- data.frame(cbind(snp.pos, phi.all))

## Plotting positon by index number
# some values that are reused on the graphing function
cex <- 0.8
num.snps <- table(snp.chrom)
labels <- c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", 
            "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", 
            "7D")[c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3))]
colours <- colors()[c(554, 144, 258, 53, 114, 450, 115)]
ylim <- c(0, 1)
## plotting
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\loci\\amova\\phi_pamk.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 150)
par(mfrow = c(3,7), oma = c(3,5,3,2), mar = c(3, 0, 0, 0))
# counts how many times tapply is called
count = 1
blah <- by(pos.phi, snp.chrom, function (x) {
  par(mar = c(2, 0, 0, 0))
  # the first graph in each row
  if (count %% 7 == 1) {
    plot(x, pch = 20, ylim = ylim, xaxt = "n", cex.axis = cex, 
         col = colours[1 + count %% 7])
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
    mtext(1, text = paste(labels[count], ": ", num.snps[count]),
                          line = 0.5, cex = cex)
  }
  # the last 6 graphs in each row
  if (count %% 7 != 1) {
    plot(x, pch = 20, ylim = ylim, xaxt = "n", yaxt = "n", cex.axis = cex, 
         col = colours[1 + count %% 7])
    abline(h = signif, col = "blue")
    abline(h = signif2, col = "red")
    mtext(1, text = paste(labels[count], ": ", num.snps[count]),
                          line = 0.5, cex = cex)
  }
  count <<- count + 1
})
title(xlab = "Marker by Chromosome Position", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5)
title("AMOVA Phi Statistic By Locus Between Pam Groups", outer = T, cex.main = 1.5)
dev.off()