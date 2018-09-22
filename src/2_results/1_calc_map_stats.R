library(tidyverse)
library(GGally)
library(SNPRelate)
library(extrafont)

source("src\\R_functions\\funcs_calc_stats.R")

gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds"
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
diff <- list(A = list(), B = list(), D = list())
num_by_segment <- by(data, snp_chrom, function (chrom) {
  if (chrom$chrom[1] %in% seq(1, 19, 3)) {
    leng$A <<- c(leng$A, max(chrom$pos))
    num$A <<- c(num$A, length(chrom$pos)) # number of markers on group
    diff$A[[chrom$chrom[1]]] <<- diff(chrom$pos)
  } else if (chrom$chrom[1] %in% seq(2, 20, 3)) {
    leng$B <<- c(leng$B, max(chrom$pos))
    num$B <<- c(num$B, length(chrom$pos))
    diff$B[[chrom$chrom[1]]] <<- diff(chrom$pos)
  } else {
    leng$D <<- c(leng$D, max(chrom$pos))
    num$D <<- c(num$D, length(chrom$pos))
    diff$D[[chrom$chrom[1]]] <<- diff(chrom$pos)
  }
})
# num SNPs and length by genome
sum(num$A)
sum(leng$A)
sum(num$B)
sum(leng$B)
sum(num$D)
sum(leng$D)
# Kength of genome covered by markers
sum(unlist(leng))
## average distance between markers overall and by genome
length(which(unlist(diff) < 0.5))
mean(unlist(diff))
mean(unlist(diff$A))
mean(unlist(diff$B))
mean(unlist(diff$D))
log10(unlist(diff))
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
lapply(diff$B, max)

# histograms and boxplots depicting the distribution of gaps on each genome
diff_log10 <- tibble(
  Genome = factor(
    c(rep("A", length(unlist(diff$A))),
      rep("B", length(unlist(diff$B))),
      rep("D", length(unlist(diff$D))),
      rep("All", length(unlist(diff)))
    ),
    levels = c("A", "B", "D", "All")
  ),
  diffs = log10(
    c(unlist(diff$A), unlist(diff$B), unlist(diff$D), unlist(diff))
  )
)
head(diff$A[[1]][2])
mean(diff_log10$diffs, na.rm = TRUE)
plots <- list()
plots[[2]] <- ggplot(diff_log10, aes(diffs, colour = Genome)) +
  geom_freqpoly() +
  scale_color_manual(values = brewer.pal(4, "Dark2")) +
  xlim(min(diff_log10$diffs), max(diff_log10$diffs))
diff_log10$Genome <- factor(diff_log10$Genome, rev(levels(diff_log10$Genome)))
plots[[1]] <- ggplot(diff_log10, aes(Genome, diffs, colour = Genome)) +
  geom_boxplot() +
  ylim(min(diff_log10$diffs), max(diff_log10$diffs)) +
  coord_flip() +
  scale_color_manual(values = rev(brewer.pal(4, "Dark2")))

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots, nrow = 2, ncol = 1, xlab = "Log10 Transformed Gap Distance",
  yAxisLabels = c("a) Boxplots", "b) Frequency Plots"),
  title = "Figure 3: Log 10 Gap Distances Pruned Map",
  legend = c(2, 1)
)

# plot the matrix
png("Results\\gaps\\gaps_dist_full.png", family = "Times New Roman",
    width = 100, height = 143, pointsize = 12, units = "mm", res = 300)
plots_matrix + theme(legend.position = "bottom")
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
    return(mean(abs(as.dist(chrom)), na.rm = T))
})

# overall avvera pairwise ld, and by genome
mean(unlist(pairwise_ld_by_chrom), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(1, 19, 3)]), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(2, 20, 3)]), na.rm = TRUE)
mean(unlist(pairwise_ld_by_chrom[seq(3, 21, 3)]), na.rm = TRUE)