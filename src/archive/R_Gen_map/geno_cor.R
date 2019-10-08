library(SNPRelate)
library(devEMF)
library(extrafont)
library(poppr)

## Setting up
# Loading the data
wheat.data <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\full_data_2.gds")
genotypes <- read.gdsn(index.gdsn(wheat.data, "genotype"))
snp.chr <- factor(read.gdsn(index.gdsn(wheat.data, "snp.chromosome")), levels = 0:21)
snp.pos <- read.gdsn(index.gdsn(wheat.data, "snp.position"))/1000
snp.id <- read.gdsn(index.gdsn(wheat.data, "snp.id"))
# Removing SNPs with high MR and low MAF
snpset <- snpgdsSelectSNP(wheat.data, autosome.only = F, maf = 0.05, missing.rate = 0.1)
snpset.id <- unlist(snpset)
informative <- match(snpset.id, snp.id)
geno.inf <- genotypes[informative,]
snp.chr.inf <- snp.chr[informative]
snp.pos.inf <- snp.pos[informative]
snpgdsClose(wheat.data)

# creating vectors without unplaced SNPs
geno.inf.non.na <- geno.inf[-which(snp.chr.inf == 0),]
snp.pos.inf.non.na <- snp.pos.inf[-which(snp.chr.inf == 0)]
snp.chr.inf.non.na <- factor(snp.chr.inf[-which(snp.chr.inf == 0)], 
                             levels = c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3)))

load("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Genotypes\\snp_indices.RData")
# Removing SNPs with high MR and low MAF
geno.kept <- genotypes[kept.indices,]
snp.chr.kept <- factor(snp.chr[kept.indices],
                       levels = c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3)))

cor.chr <- by(geno.inf.non.na, snp.chr.inf.non.na, function (x) {
  cor(t(x))
})


# Chromosome labels
labels.reordered <- c("NA", "1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D",
                      "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", 
                      "7A", "7B", "7D")[c(seq(2, 20, 3), seq(3, 21, 3), seq(4, 22, 3), 1)]

#####################
## informative NON NA SNPs
num.snps <- table(snp.chr.inf.non.na)
png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\cor\\cor_all.png", 
    family="Times New Roman", width = 23, height = 12, pointsize = 20, units = "in", res = 300)
par(mfrow = c(3,7), oma = c(3,5,3,10), mar = c(3, 0, 0, 0))
count <- 1
cex = 0.8
foo <- tapply(snp.pos.inf.non.na, snp.chr.inf.non.na, function (x) {
#   set.seed(3000)
#   x <- sort(jitter(x, amount = 0.1))
  if (count %% 7 == 1) {
    image(abs(cor.chr[[count]]), col = topo.colors(4), xaxt = "n")
    mtext(1, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = 1, cex = cex)
  }
  if (count %% 7 != 1) {
    image(abs(cor.chr[[count]]), col = topo.colors(4), xaxt = "n", yaxt = "n")
    mtext(1, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = 1, cex = cex)
  }
  count <<- count + 1
})
# titles and axes labels
title("Fig 9: Heatmaps of Absolute Correlation between Markers on each Chromosome", out = T, cex = 1.5)
title(xlab = "Markers Ordered by Genetic Map", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "Markers Ordered by Genetic Map", outer = T, cex.lab = 1.5)
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
add_legend(0.95, 0.99, legend=c("0.00-0.24", "0.25-0.49", "0.50-0.74", "0.75-1.00"),
           pch=15, col=topo.colors(4), cex = 1.2)
dev.off()


######################
## pruned
cor.chr <- by(geno.kept, snp.chr.kept, function (x) {
  cor(t(x))
})

num.snps <- table(snp.chr.kept)

png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Figures\\cor\\cor_kept.png", 
    family="Times New Roman", width = 12, height = 6, pointsize = 10, units = "in", res = 150)
par(mfrow = c(3,7), oma = c(3,5,3,10), mar = c(3, 0, 0, 0))
count <- 1
for(chr in names(cor.chr)) {
  if (count %% 7 == 1) {
    image(abs(cor.chr[[chr]]), col = topo.colors(4), xaxt = "n")
    mtext(1, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = 1, cex = cex)
  }
  if (count %% 7 != 1) {
    image(abs(cor.chr[[chr]]), col = topo.colors(4), xaxt = "n", yaxt = "n")
    mtext(1, text = paste(labels.reordered[count], ": ", num.snps[count]),
          line = 1, cex = cex)
  }
  count <<- count + 1
}
# titles and axes labels
title("Heatmaps of Absolute Correlation between Markers on each Chromosome", out = T)
title(xlab = "Markers", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "Markers", outer = T, cex.lab = 1.5)
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
add_legend(0.92, 0.95, legend=c("0.00-0.24", "0.25-0.49", "0.50-0.74", "0.75-1.00"),
       pch=15, col=topo.colors(4), cex = 1.2)
dev.off()