library(extrafont)
library(SNPRelate)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\R\\functions.R")

## useful data
load("Data\\Intermediate\\Aligned_genes\\wheat_genes_contigs.RData")
wheat <-
  snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
snp.chrom <- read.gdsn(index.gdsn(wheat, "snp.chromosome"))
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
sample.id <- read.gdsn(index.gdsn(wheat, "sample.id"))

# get the pop code
desig <-
  read.gdsn(index.gdsn(wheat, "samp.annot/designation"))
year <-
  read.gdsn(index.gdsn(wheat, "samp.annot/Year"))

## for making subsets of the samples for the different categories excluding unknown samples
index.hrs <-
  c(which(desig == "HRS"))
index.hrs_hrw <-
  c(index.hrs, which(desig == "HRW"))
young_vs_old <- vector(length = 234)
young_vs_old[which(as.integer(as.character(year[index.hrs])) <= 1986)] <- "old"
young_vs_old[which(as.integer(as.character(year[index.hrs])) > 1986)] <- "young"

# fst.all <- snpgdsFst(wheat, sample.id = sample.id[index.hrs_hrw], population = factor(desig[index.hrs_hrw]), method = "W&H02", autosome.only = F)$FstSNP
fst.all <- snpgdsFst(wheat, sample.id = sample.id[index.hrs], population = factor(young_vs_old), method = "W&H02", autosome.only = F, maf = NaN, missing.rate = NaN, remove.monosnp = F)$FstSNP
snpgdsClose(wheat)

## finding phi results for each locus and the top quantiles
signif <- quantile(fst.all, prob = 0.995, na.rm = T)
mean.fst.all.list <-
  by(cbind(snp.pos, fst.all), snp.chrom, function (x) {
    sw_calc(x[, 1], x[, 2])
  })
mean.fst.all <- vector()
for (i in 1:length(mean.fst.all.list)) {
  mean.fst.all <-
    c(mean.fst.all, mean.fst.all.list[[i]][-which(mean.fst.all.list[[i]][, 2] == 0), 2])
}
signif2 <- quantile(mean.fst.all, prob = 0.995, na.rm = T)

## plotting
png(
  "Results\\loci\\full\\full_fst_hrs_vs_hrw.png",
  family = "Times New Roman",
  width = 6,
  height = 6,
  pointsize = 10,
  units = "in",
  res = 300
)
plots(snp.pos, fst.all, genes_contigs, snp.chrom)
title(
  xlab = "Marker by Chromosome Position",
  outer = T,
  cex.lab = 1.5,
  line = 0.5
)
title(ylab = "AMOVA Phi Statistic",
      outer = T,
      cex.lab = 1.5)
title(
  "AMOVA Phi Statistic By Locus Between HRS and HRW Wheats",
  outer = T,
  cex.main = 1.5
)
dev.off()
