library(SNPRelate)
library(extrafont)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("src\\analysis\\R\\R_functions\\funcs_draw_loci.R")

## genes split up by contig
mainGenes <- get(load("Data\\Intermediate\\Aligned_genes\\top_main_genes_contigs.RData"))
resiGenes <- get(load("Data\\Intermediate\\Aligned_genes\\top_resistance_genes_contigs.RData"))
mainGenes$Pos <- mainGenes$Pos/1000000
resiGenes$Pos <- resiGenes$Pos/1000000

# prints out the eh values of the extreme phi markers of each group in the compariosn on the smae graph with the second group as negative values
for (subset in list(c("HRS", "SWS", "HRS_SWS"), c("HRS", "HRW", "HRS_HRW"), c("HRW", "SWS", "HRW_SWS"))) {
  gdsSubset <- paste0("Data\\Intermediate\\GDS\\", subset[1], "_phys_subset_sample.gds")
  source("src\\R_functions\\data_loading.R")
  group1EH <- EH(genotypes)
  gdsSubset <- paste0("Data\\Intermediate\\GDS\\", subset[2], "_phys_subset_sample.gds")
  source("src\\R_functions\\data_loading.R")
  group2EH <- -EH(genotypes)
  snp.pos <- snp.pos/1000000
  
  EHs <- cbind(group1EH, group2EH)
  
  amova <- get(load(paste0("Data\\Intermediate\\Adegenet\\amova_", subset[3], ".RData")))
  phiAll <- phi.all(amova)
  signif <- quantile(phiAll, prob = 0.975, na.rm = T)
  
  EHs[which(phiAll < signif),] <- NA
  
  # plotting
  png(paste0("Results\\loci\\EH\\EH_between_", subset[1], "_and_", subset[2], ".png"),
      family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
  plots2(snp.pos, EHs, mainGenes, resiGenes, snp.chrom, ylim = c(-0.5, 0.5))
  title(xlab = "Marker Position in Mb", outer = T, cex.lab = 1.5, line = 2.5)
  title(ylab = "Expectged Heterozygosity", outer = T, cex.lab = 1.5, line = 0)
  if (subset[1] == "HRS") subset[1] <- "CHRS"
 
  par(xpd = NA)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  title(main = paste("Difference in Expected Heterozygosity Between", subset[1], "and", subset[2], "At Markers With Extreme Phi Values"),
        outer = T, cex.main = 1.5, line = -1.5)
  legend(0.39, -1.01, legend = c("Designation and Dwarfing Genes", "Resistance Genes"),
         pch = 19, col = c(colourSet[15], colourSet[19]), cex = 0.6, horiz = T)
  dev.off()
}

## calculates difference between eh values of comparisons and plots them (loses magnitude)
# load data
# for (subset in list(c("HRS", "SWS", "HRS_SWS"), c("HRS", "HRW", "HRS_HRW"), c("HRW", "SWS", "HRW_SWS"))) {
#   gdsSubset <- paste0("Data\\Intermediate\\GDS\\", subset[1], "_phys_subset_sample.gds")
#   source("Analysis\\R\\R_functions\\data_loading.R")
#   group1EH <- EH(genotypes)
#   gdsSubset <- paste0("Data\\Intermediate\\GDS\\", subset[2], "_phys_subset_sample.gds")
#   source("Analysis\\R\\R_functions\\data_loading.R")
#   group2EH <- EH(genotypes)
#   snp.pos <- snp.pos/1000000
# 
#   EHDiff <- group1EH - group2EH
# 
#   amova <- get(load(paste0("Data\\Intermediate\\Adegenet\\amova_", subset[3], ".RData")))
#   phiAll <- phi.all(amova)
#   signif <- quantile(phiAll, prob = 0.975, na.rm = T)
# 
#   EHDiff[which(phiAll < signif)] <- NA
# 
#   # plotting
#   png(paste0("Results\\loci\\EH\\EH_between_", subset[1], "_and_", subset[2], ".png"),
#       family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
#   plots2(snp.pos, EHDiff, mainGenes, resiGenes, snp.chrom, ylim = c(-0.5, 0.5))
#   title(xlab = "Marker Position in Base Pairs", outer = T, cex.lab = 1.5, line = 2.5)
#   title(ylab = "Expectged Heterozygosity", outer = T, cex.lab = 1.5, line = 0)
#   title(main = paste("Difference in Expected Heterozygosity Between", subset[1], "and", subset[2], "At Markers With Extreme Phi Values"), outer = T, cex.main = 1.25, line = 1)
# 
#   par(xpd = NA)
#   par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
#   plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#   legend(0.39, -1.01, legend = c("Designation and Dwarfing Genes", "Resistance Genes"),
#          pch = 19, col = c(colourSet[15], colourSet[19]), cex = 0.6, horiz = T)
#   dev.off()
# }

gdsSubset <- paste0("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
source("src\\analysis\\R\\R_functions\\data_loading.R")
snp.pos <- snp.pos/1000000

haplotypes <- as.list(rep(NA, 21))

A1 <- which(snp.chrom == "1")
haplotypes[[1]] <- match(which(snp.pos[A1] > 7E7), which(snp.pos[A1] < 3.2E8))
haplotypes[[1]] <- haplotypes[[1]][!is.na(haplotypes[[1]])]

A2 <- which(snp.chrom == "4")
haplotypes[[4]] <- match(which(snp.pos[A2] > 2.1E8), which(snp.pos[A2] < 4.7E8))
haplotypes[[4]]<- haplotypes[[4]][!is.na(haplotypes[[4]])]

A4 <- which(snp.chrom == "10")
haplotypes[[10]] <- match(which(snp.pos[A4] > 2.3E8), which(snp.pos[A4] < 4.6E8))
haplotypes[[10]] <- haplotypes[[10]][!is.na(haplotypes[[10]])]

B5 <- which(snp.chrom == "14")
haplotypes[[14]] <- match(which(snp.pos[B5] > 1.1E8), which(snp.pos[B5] < 2.1E8))
haplotypes[[14]] <- haplotypes[[14]][!is.na(haplotypes[[14]])]

A6 <- which(snp.chrom == "16")
haplotypes[[16]] <- match(which(snp.pos[A6] > 1.7E8), which(snp.pos[A6] < 4.45E8))
haplotypes[[16]] <- haplotypes[[16]][!is.na(haplotypes[[16]])]

B6 <- which(snp.chrom == "17")
haplotypes[[17]] <- match(which(snp.pos[B6] > 2.5E8), which(snp.pos[B6] < 3.8E8))
haplotypes[[17]] <- haplotypes[[17]][!is.na(haplotypes[[17]])]

A7 <- which(snp.chrom == "19")
haplotypes[[19]] <- match(which(snp.pos[A7] > 3.1E8), which(snp.pos[A7] < 4.45E8))
haplotypes[[19]] <- haplotypes[[19]][!is.na(haplotypes[[19]])]

# source("Analysis\\R\\R_functions\\funcs_draw_loci.R")

## plots the eh values for all markers
for (subset in c("wheat", "HRS", "HRW", "SWS")[1]) {
  gdsSubset <- paste0("Data\\Intermediate\\GDS\\", subset, "_phys_subset_sample.gds")
  source("src\\R_functions\\data_loading.R")
  snp.pos <- snp.pos/1000000
  print(subset)
  
  # plotting
  png(paste0("Results\\loci\\EH\\", subset, "_phys_EH.png"),
      family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
  plots3(snp.pos, EH(genotypes), mainGenes, resiGenes, snp.chrom, ylim = c(0, 0.5),
         oma = c(4,0,3,4), colour = TRUE, plotGenes = FALSE, haplotypes = haplotypes)
  title(xlab = "Marker Position in Mb", outer = T, cex.lab = 1.5, line = 2.5)
  mtext("Expectged Heterozygosity", outer = T, cex.lab = 1.5, line = 2.5, side = 4)
  if (subset == "HRS") subset <- "CHRS"
  
  par(xpd = NA)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  if (subset == "wheat") {
    title(paste("Expected Heterozygosity Within All Varieties Across Chromosomes Plotted by Physical Map Position"),
          outer = T, cex.main = 1.5, line = -1.5)
  } else {
    title(paste("Expected Heterozygosity Within", subset, "Designates Across Chromosomes Plotted by Physical Map Position"),
          outer = T, cex.main = 1.5, line = -1.5)
  }
  dev.off()
}

# source("Analysis\\R\\R_functions\\funcs_draw_loci.R")

## plots the eh values for all markers
for (subset in c("wheat", "HRS", "HRW", "SWS")[1]) {
  gdsSubset <- paste0("Data\\Intermediate\\GDS\\", subset, "_gen_subset_sample.gds")
  source("src\\R_functions\\data_loading.R")
  snp.pos <- snp.pos/1000000
  print(subset)
  
  # plotting
  png(paste0("Results\\loci\\EH\\", subset, "_gen_EH.png"),
      family="Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
  plots3(snp.pos, EH(genotypes), mainGenes, resiGenes, snp.chrom, ylim = c(0, 0.5), oma = c(4,0,3,4),
         colour = TRUE, plotGenes = FALSE, haplotypes = haplotypes)
  title(xlab = "Marker Position in cM", outer = T, cex.lab = 1.5, line = 2.5)
  mtext("Expectged Heterozygosity", outer = T, cex.lab = 1.5, line = 2.5, side = 4)
  if (subset == "HRS") subset <- "CHRS"
  
  par(xpd = NA)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  if (subset == "wheat") {
    title(paste("Expected Heterozygosity Within All Varieties Across Chromosomes Plotted By Genetic Map Position"),
          outer = T, cex.main = 1.5, line = -1.5)
  } else {
    title(paste("Expected Heterozygosity Within", subset, "Designates Across Chromosomes Plotted By Genetic Map Position"),
          outer = T, cex.main = 1.5, line = -1.5)
  }
  dev.off()
}

