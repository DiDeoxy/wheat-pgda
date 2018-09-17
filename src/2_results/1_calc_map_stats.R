library(tidyverse)
library(SNPRelate)
library(extrafont)

source("src\\R_functions\\funcs_calc_stats.R")

gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")

# how many snps?
length(snp_id)

# count how often for each marker the different homozygotes and missing data
# occurs
counts <- apply(genotypes, 1, function (marker) {
  count <- as.data.frame(table(marker))
  if (nrow(count) == 3) {
    return(count[, 2])
  } else {
    return(c(count[, 2], 0))
  }
})
## find the minor allel frequency and missing rate of each marker
stats <- as.data.frame(t(apply(counts, 2, function (marker) {
                            maf <- min(marker[1:2]) / sum(marker[1:2])
                            mr <- marker[3] / sum(marker)
                            return(c(maf = maf, mr = mr))
                         })))
## mean maf and mr
mean(stats$maf)
mean(stats$mr)

## number of markers and mean distances between the genome
leng <- list(A = vector(), B = vector(), D = vector())
num <- list(A = vector(), B = vector(), D = vector())
diff <- list(A = vector(), B = vector(), D = vector())
num_by_segment <- by(data, snp_chrom, function (chrom_pos) {
  if (chrom_pos$chrom[1] %in% seq(1, 19, 3)) {
    leng$A <<- c(leng$A, max(chrom_pos$pos))
    num$A <<- c(num$A, length(chrom_pos$pos)) # number of markers on group
    diff$A <<- c(diff$A, diff(chrom_pos$pos))
  } else if (chrom_pos$chrom[1] %in% seq(2, 20, 3)) {
    leng$B <<- c(leng$B, max(chrom_pos$pos))
    num$B <<- c(num$B, length(chrom_pos$pos))
    diff$B <<- c(diff$B, diff(chrom_pos$pos))
  } else {
    leng$D <<- c(leng$D, max(chrom_pos$pos))
    num$D <<- c(num$D, length(chrom_pos$pos))
    diff$D <<- c(diff$D, diff(chrom_pos$pos))
  }
  count <<- count + 1
})
# mac position on each chromosome and genome
leng
sum(leng$A)
sum(leng$B)
sum(leng$D)
# num SNPs by chromosome and genopme
num
sum(num$A)
sum(num$B)
sum(num$D)
## average distance between markers overall and by genome
mean(unlist(diff))
mean(diff$A)
mean(diff$B)
mean(diff$D)

# find the min length of the top percentile of gaps
top_percentile <- quantile(unlist(diff), prob = 0.99, na.rm = T)
top_percentile

# the number of top percentile gaps on each genome
length(which(unlist(diff$A) >= top_percentile))
length(which(unlist(diff$B) >= top_percentile))
length(which(unlist(diff$D) >= top_percentile))
# the longest gap on each genome
max(unlist(diff$A))
max(unlist(diff$B))
max(unlist(diff$D))

# histograms and boxplots depicting the distribution of gaps on each genome
png("Results\\gaps\\gaps_dist.png", family="Times New Roman",
    width = 12, height = 6, pointsize = 20, units = "in", res = 500)
par(mar = c(0, 1, 1, 1), oma = c(4, 3, 3, 0))
nf <- layout(mat = matrix(c(1, 3, 5, 7, 2, 4, 6, 8), 2, 4, byrow = TRUE), 
             height = c(4, 1))
layout.show(nf)
for (chrom in names(diff)) {
  hist(log10(unlist(diff[[chrom]])), col = "pink", main = chrom, xlab = "",
       ylab = "", xaxt = "n")
  boxplot(log10(unlist(diff[[chrom]])), horizontal = TRUE, outline = TRUE, 
          frame = F, col = "green1", width = 10)
}
hist(log10(unlist(diff)), col = "pink", main = "All", xlab = "", ylab = "",
     xaxt = "n")
boxplot(log10(unlist(diff)), horizontal = TRUE, outline = TRUE,
        frame = F, col = "green1", width = 10)
title(xlab = "Log 10 Power", outer = T, cex.lab = 1.5, line = 2.5)
title(ylab = "Count", outer = T, cex.lab = 1.5, line = 1.5)
title(paste("Distributions of Distances Between Markers by Genome"), outer = T,
      cex.main = 1.5, line = 1)
dev.off()

# Calcualte ld between all markers on each chromosome
full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_subset_sample.gds")
ld_mats <- by(data, snp_chrom, function (chrom) {
  abs(snpgdsLDMat(full, method = "composite", snp.id = chrom$id,
                  slide = -1)$LD) 
})
snpgdsClose(full)

# calculate the average pairwise ld by chromosome
pairwise_ld_by_chrom <- lapply(ld_mats, function (chrom) {
    return(mean(as.dist(data$LD), na.rm = T))
    # return(as.vector(as.dist(chrom$LD^2)))
})

# overall avvera pairwise ld, and by genome
mean(unlist(pairwise_ld_by_chrom), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(1, 19, 3)]), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(2, 20, 3)]), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(3, 21, 3)]), na.rm = TRUE)