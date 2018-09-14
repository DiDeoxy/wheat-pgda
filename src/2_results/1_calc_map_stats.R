library(tidyverse)
library(SNPRelate)
library(extrafont)

source("src\\R_functions\\funcs_calc_stats.R")

gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")

mb <- 1000000

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
num <- data.frame(A = 0, B = 0, D = 0)
diff <- list(A = list(), B = list(), D = list())
count <- 1
num_by_segment <- by(snp_pos, snp_chrom, function (pos) {
  if (count %in% seq(1, 19, 3)) {
    num$A <<- sum(num$A, length(pos)) # number of markers on group
    if (count == 1) {
      diff$A[[1]] <<- diff(pos) # distances between markers
    } else {
      diff$A[[(count %/% 3) + 1]] <<- diff(pos) # for each chromosome by genome
    }
  } else if (count %in% seq(2, 20, 3)) {
    num$B <<- sum(num$B, length(pos))
    if (count == 1) {
      diff$B[[1]] <<- diff(pos)
    } else {
      diff$B[[(count %/% 3) + 1]] <<- diff(pos)
    }
  } else {
    num$D <<- sum(num$D, length(pos))
    diff$D[[count / 3]] <<- diff(pos)
  }
  count <<- count + 1
})
# num SNPs by genome
num
## average distance between markers overall and by genome
mean(unlist(diff)) / mb
mean(unlist(diff$A)) / mb
mean(unlist(diff$B)) / mb
mean(unlist(diff$D)) / mb

# histograms and boxplots depicting the distribution of gaps on each genome
png("Results\\gaps\\gaps_dist.png", family="Times New Roman",
    width = 12, height = 6, pointsize = 20, units = "in", res = 500)
par(mar = c(0, 1, 1, 1), oma = c(4, 3, 3, 0))
nf <- layout(mat = matrix(c(1, 3, 5, 7, 2, 4, 6, 8), 2, 4, byrow = TRUE), 
             height = c(4,1))
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

# find the min length of the top percentile of gaps
top_percentile <- quantile(unlist(diff), prob = 0.99, na.rm = T)
top_percentile

# the number of top percentile gaps on each genome
length(which(unlist(diff$A) >= top_percentile))
length(which(unlist(diff$B) >= top_percentile))
length(which(unlist(diff$D) >= top_percentile))
# the longest gap on each genome
max(unlist(diff$A)) / mb
max(unlist(diff$B)) / mb
max(unlist(diff$D)) / mb

# Calcualte ld between all markers on each chromosome
wheat <- snpgdsOpen(
    "Data\\intermediate\\GDS\\full_phys_subset_sample.gds")
chroms <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste,
                            sep = "")))
ld_mat <- list()
for (i in unique(snp_chrom)) {
    ld_mat[[chroms[i]]] <- c(snpgdsLDMat(wheat, method = "composite",
                                snp.id = snp_id[which(snp_chrom == i)],
                                slide = -1),
                             list(pos = snp_pos[which(snp_chrom == i)]))
}
snpgdsClose(wheat)

# calculate the average pairwise ld by chromosome
pairwise_ld_by_chrom <- lapply(ld_mat, function (chrom) {
    # return(mean(abs(as.dist(data$LD)), na.rm = T))
    return(as.vector(as.dist(chrom$LD^2)))
})

# overall avvera pairwise ld, and by genome
mean(unlist(pairwise_ld_by_chrom), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(1, 19, 3)]), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(2, 20, 3)]), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(3, 21, 3)]), na.rm = TRUE)