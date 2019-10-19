library(SNPRelate)
library(extrafont)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
load("Data\\Intermediate\\Aligned_genes\\HRS_genes_contigs.RData")
source("Analysis\\R\\functions.R")

## useful data
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\HRS_phys_subset_sample.gds")
years <- read.gdsn(index.gdsn(wheat, "samp.annot/Year"))
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chrom <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
summary(years)
fst_years <- list()
for (year in seq(1980, 2000, 1)) {
  samp.comp <- integer(length(years))
  samp.comp[which(as.numeric(as.character(years)) <= year)] <- 1
  
  fst_years[[as.character(year)]] <- snpgdsSlidingWindow(wheat, FUN = "snpgdsFst", unit = "locus", winsize = 1, shift = 1, 
                                      population = factor(samp.comp), sample.id = sample.id, 
                                      method ="W&H02", win.start = 0, as.is = "numeric", remove.monosnp = F, maf = NaN, MR = NaN)
}
snpgdsClose(wheat)


## plotting
for (year in as.character(seq(1980, 2000, 1))) {
  ## finding fst quantiles
  fst_year <- vector()
  for (i in c(seq(3, 167, 4))) {
    fst_year <- c(fst_year, fst_years[[year]][i])
  }
  fst_year <- unlist(fst_year)
  signif <- quantile(fst_year, prob = 0.95, na.rm = T)
  signif2 <- quantile(fst_year, prob = 0.99, na.rm = T)
  
  png(paste("Results\\loci\\hrs_amova_pre_vs_post_", year, ".png", sep = ""),
      family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
  plots(snp.pos, fst_year, genes_contigs, snp.chrom, ylim = c(0, 0.5))
  title(xlab = "Marker by Chromosome Position", outer = T, cex.lab = 1.5, line = 0.5)
  title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5)
  title(paste("AMOVA Phi Statistic By Locus Between pre vs post", year), outer = T, cex.main = 1.5)
  par(xpd=NA)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend=c("Raw 0.5% Threshold",
                            "Smoothed 0.5% Threshold"),
         pch=19, col=c("red", "blue"), cex = 0.6, horiz = T)
  dev.off()
}