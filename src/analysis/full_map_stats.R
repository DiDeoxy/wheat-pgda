library(SNPRelate)
library(extrafont)
# install.packages("GeneticSubsetter")
# library(GeneticSubsetter)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\R\\functions\\funcs_calc_stats.R")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

# load("Data\\Intermediate\\Aligned_genes\\top_main_genes_contigs.RData")
dim(genotypes)
##
# genotypes <- replace(genotypes, genotypes == 0, -1)
# genotypes <- replace(genotypes, genotypes == 2, 1)
# genotypes <- replace(genotypes, genotypes == 3, 0)

# PICs <- vector()
# for (i in 1:dim(genotypes)[1]) {
#   PICs <- c(PICs, PicCalc(t(genotypes[i,])))
# }
# signif <- quantile(PICs, prob = 0.95, na.rm = T)
# signif2 <- quantile(PICs, prob = 0.99, na.rm = T)
# mean(PICs)

# png("Analysis\\Figures\\loci\\pic_hrs.png", 
#     family="Calibri", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
# plots2(snp.pos, PICs, genes_contigs, snp.chrom)
# dev.off()

counts <- apply(genotypes, 1, function (x) {
  count <- as.data.frame(table(x))
  if (dim(count)[1] == 3) {
    return(count[,2])
  } else {
    return(c(count[,2], 0))
  }
})
stats <- apply(counts, 2, function (x) {
  maf <- min(x[1:2])/365
  mr <- x[3]/365
  return(c(maf = maf, mr = mr))
})
## MAF and MR
mean(stats[1,])
mean(stats[2,])

## numbers and mean distances
num <- data.frame(A = 0, B = 0, D = 0)
diff <- list(A = list(), B = list(), D = list())
count <- 1
i <- 1
num_by_segment <- by(snp.pos, snp.chrom, function (x) {
  if (count %in% seq(1, 19, 3)) {
    num$A <<- sum(num$A, length(x))
    if (count == 1) {
      diff$A[[1]] <<- diff(x)
    } else {
      diff$A[[(count %/% 3) + 1]] <<- diff(x)
    }
    
  } else if (count %in% seq(2, 20, 3)) { 
    num$B <<- sum(num$B, length(x))
    if (count == 1) {
      diff$B[[1]] <<- diff(x)
    } else {
      diff$B[[(count %/% 3) + 1]] <<- diff(x)
    }
  } else {
    num$D <<- sum(num$D, length(x))
    diff$D[[count / 3]] <<- diff(x)
  }
  count <<- count + 1
})
## num SNPs
num

## average dist
mean(unlist(diff))
mean(unlist(diff$A))
mean(unlist(diff$B))
mean(unlist(diff$D))

png("Results\\gaps_dist_ldpruned.png", family="Times New Roman", width = 12, height = 6, pointsize = 20, units = "in", res = 500)
par(mar=c(0, 1, 1, 1), oma = c(4, 3, 3, 0))
nf <- layout(mat = matrix(c(1,3,5,7,2,4,6,8), 2, 4, byrow=TRUE), height = c(4,1))
layout.show(nf)
for (chrom in names(diff)) {
  hist(log10(unlist(diff[[chrom]])), col = "pink", main = chrom, xlab = "", ylab = "", xaxt = "n")
  boxplot(log10(unlist(diff[[chrom]])), horizontal=TRUE, outline=TRUE,
          frame=F, col = "green1", width = 10)
}
hist(log10(unlist(diff)), col = "pink", main = "All", xlab = "", ylab = "", xaxt = "n")
boxplot(log10(unlist(diff)), horizontal=TRUE, outline=TRUE,
        frame=F, col = "green1", width = 10)
title(xlab = "Log 10 Power", outer = T, cex.lab = 1.5, line = 2.5)
title(ylab = "Count", outer = T, cex.lab = 1.5, line = 1.5)
title(paste("Distributions of Distances Between Markers by Genome"), outer = T, cex.main = 1.5, line = 1)
dev.off()

## Gaps

longest <- quantile(unlist(diff), prob = 0.99, na.rm = T)
longest
length(which(unlist(diff$A) >= longest))
length(which(unlist(diff$B) >= longest))
length(which(unlist(diff$D) >= longest))
order(unlist(diff$A), decreasing = T)
unlist(diff$A)[2199]
max(unlist(diff$B))
max(unlist(diff$D))

## LD
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
# labels <- as.vector(t(outer(chromes, c("part1", "part2"), paste, sep="_")))

wheat <- snpgdsOpen("Data\\intermediate\\GDS\\wheat_phys_subset_both.gds")
ld_mat <- list()
for(i in levels(as.factor(snp.chrom))) {
  ld_mat[[labels[as.numeric(i)]]] <- c(snpgdsLDMat(wheat, method = "composite", 
                                                   snp.id = snp.id[which(snp.chrom == i)], 
                                                   slide = -1),
                                       list(pos = snp.pos[which(snp.chrom == i)]))
}
snpgdsClose(wheat)


pairwise_ld <- lapply(ld_mat, function (data) {
  ld <- vector()
  # ld <- c(ld, mean(abs(as.dist(data$LD)), na.rm = T))
  ld <- c(ld, mean(as.dist(data$LD^2), na.rm = T))
  return(ld)
})

genome_ld <- list(A = list(), B = list(), D = list())
count <- 1
blah <- lapply(pairwise_ld, function (pairs) {
  if (count %in% seq(1, 19, 3)) {
    if (count == 1) {
      genome_ld$A[[1]] <<- pairs
    } else {
      genome_ld$A[[(count %/% 3) + 1]] <<- pairs
    }
  } else if (count %in% seq(2, 20, 3)) { 
    if (count == 1) {
      genome_ld$B[[1]] <<- pairs
    } else {
      genome_ld$B[[(count %/% 3) + 1]] <<- pairs
    }
  } else {
    genome_ld$D[[count / 3]] <<- pairs
  }
  count <<- count + 1
})
mean(unlist(genome_ld))

lapply(genome_ld, function (x) {
  mean(unlist(x))
})


mean(unlist(pairwise_ld), na.rm = T)

