library(SNPRelate)
library(extrafont)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\R\\functions\\funcs_draw_loci.R")
source("Analysis\\R\\functions\\HRS_subset_sample_data_loading.R")
source("Analysis\\R\\functions\\colour_sets.R")
source("Analysis\\R\\functions\\funcs_calc_stats.R")

## loading pre-calcualted per locus amova results
load("Data\\Formatted\\amova_psr_hrs.RData")
load("Data\\HRS_genes_contigs.RData")

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
png("Analysis\\Figures\\loci\\HRS_amova_Group1_vs_Group4.png", 
    family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
plots(snp.pos, phi.all, genes_contigs, snp.chrom)
title(xlab = "Marker by Chromosome Position", outer = T, cex.lab = 1.5, line = 0.5)
title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5)
title("AMOVA Phi Statistic By Locus Between Group 1 and Group 4 Wheats", outer = T, cex.main = 1.5)
dev.off()
