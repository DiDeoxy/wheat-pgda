library(extrafont)
library(SNPRelate)

setwd("C:\\Users\\Max_H.DESKTOP-AJ57KB6\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\Scripts\\Physical map\\functions.R")

## useful data
wheat <- snpgdsOpen("Data\\Formatted\\wheat_phys_subset_sample.gds")
snp.chrom <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snpgdsClose(wheat)
## loading pre-calcualted per locus amova results
load("Data\\Formatted\\amova_hrs_hrw.RData")
phi.all1 <- phi.all(wheat.amova)
load("Data\\Formatted\\amova_hrs_sws.RData")
phi.all2 <- phi.all(wheat.amova)

load("Data\\wheat_genes_contigs.RData")

phi.all <- phi.all2 - phi.all1

## plotting
png("Analysis\\Figures\\loci\\full_amova_hrs_dif.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
plots3(snp.pos, phi.all, genes_contigs, snp.chrom)
title(xlab = "Marker by Chromosome Position", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5)
title("AMOVA Phi Statistic By Locus Difference Between HRS vs HRW and HRS vs SWS Wheats", outer = T, cex.main = 1.5)
dev.off()
