library(SNPRelate)
library(extrafont)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\Scripts\\Physical map\\functions.R")

## useful data
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\HRS_phys_subset_sample.gds")
snp.chrom <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snpgdsClose(wheat)
## loading pre-calcualted per locus amova results
load("Data\\Intermediate\\Adegenet\\amova_hrs2_hrs3.RData")
load("Data\\Aligned_genes\\HRS_genes_contigs.RData")

## finding phi results for each locus and the top quantiles
phi.all <- phi.all(wheat.amova)
signif <- quantile(phi.all, prob = 0.995, na.rm = T)
mean.phi.all.list <- by(cbind(snp.pos, phi.all), snp.chrom, function (x) { swpos(x[,1], x[,2]) })
mean.phi.all <- vector()
for (i in 1:length(mean.phi.all.list)) {
  mean.phi.all <- c(mean.phi.all, mean.phi.all.list[[i]][-which(mean.phi.all.list[[i]][,2] == 0),2])
}
signif2 <- quantile(mean.phi.all, prob = 0.995, na.rm = T)

## plotting
png("Analysis\\Figures\\loci\\HRS_amova_Group2_vs_Group3.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
plots(snp.pos, phi.all, genes_contigs, snp.chrom)
length(snp.chrom)
title(xlab = "Marker by Chromosome Position", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5)
title("AMOVA Phi Statistic By Locus Between Group 2 and Group 3 Wheats", outer = T, cex.main = 1.5)
dev.off()
