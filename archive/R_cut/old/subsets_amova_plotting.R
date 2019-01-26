library(extrafont)
library(SNPRelate)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA")
source("Analysis\\R\\functions\\funcs_draw_loci.R")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
source("Analysis\\R\\functions\\data_loading.R")
which(snp.id == "BobWhite_c1740_97")
snp.chrom[5654]
snp.pos[5654]

# aligned genes
mainGenes <- get(load("Data\\Intermediate\\Aligned_genes\\top_main_genes_contigs.RData"))
resiGenes <- get(load("Data\\Intermediate\\Aligned_genes\\top_resistance_genes_contigs.RData"))
knownGenes <- get(load("Data\\Intermediate\\Aligned_genes\\known_genes_groups.RData"))

for (group in c("CCP", "HRS_SWS", "HRS_HRW", "Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")[1]) {

  amova <- get(load(paste0("Data\\Intermediate\\Adegenet\\amova_", group, ".RData")))
  phiAll <- phi.all(amova)
  signif <- quantile(phiAll, prob = 0.995, na.rm = T)
  mean.phi.all.list <-
    by(cbind(snp.pos, phiAll), snp.chrom, function (x) {
      sw_calc(x[, 1], x[, 2])
    })
  mean.phi.all <- vector()
  for (i in 1:length(mean.phi.all.list)) {
    mean.phi.all <-
      c(mean.phi.all, mean.phi.all.list[[i]][-which(mean.phi.all.list[[i]][, 2] == 0), 2])
  }
  signif2 <- quantile(mean.phi.all, prob = 0.995, na.rm = T)

  ## plotting
  png(paste0("Results\\loci\\amova\\full_amova_", group, ".png"),
      family = "Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
  plots(snp.pos, phiAll, mainGenes, resiGenes, snp.chrom, raw = FALSE, colour = TRUE, ab = "smooth")
  title(xlab = "Marker Position in Base Pairs", outer = T, cex.lab = 1.5, line = 2.5)
  title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5)
  names <- unlist(strsplit(group, "_"))
  if (length(names) == 2) {
    title(paste("AMOVA Phi Statistic By Locus Between", names[1], "and", names[2], "Wheats"), outer = T, cex.main = 1.5)
  } else {
    title(paste("AMOVA Phi Statistic By Locus Between", group, "and", tolower(group), "Wheats"), outer = T, cex.main = 1.5)
  }

  par(xpd = NA)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  # legend(0.39, -1.01, legend = c("Raw 0.5% Threshold", "Smoothed 0.5% Threshold"), pch = 19, col = c("red", "blue"), cex = 0.6, horiz = T)
  legend(0.37, -1.01, legend = c("Top 0.5% Threshold"), pch = 19, col = c("blue"), cex = 0.6, horiz = T)
  legend(-0.97, -1.01, legend = c("Designation and Dwarfing Genes", "Resistance Genes"), 
         pch = 19, col = c(colourSet[15], colourSet[19]), cex = 0.6, horiz = T)
  dev.off()
}
# 
# png(paste0("Results\\loci\\amova\\full_amova_genes.png"), family = "Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
# par(mfrow = c(5,1), mar = c(1,4,1,2), oma = c(4,0,0,0))
# emptyFrame <- data.frame(one = rep(NA, 15115), two = rep(NA, 15115), three = rep(NA, 15115))
# for (group in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
#   amova <- get(load(paste0("Data\\Intermediate\\Adegenet\\amova_", group, ".RData")))
#   phiAll <- phi.all(amova)
#   signif <- quantile(phiAll, prob = 0.999, na.rm = T)
# 
#   ## plotting
#   count <- 1
#   by(cbind(snp.pos, phiAll, emptyFrame, knownGenes), snp.chrom, function (chr) {
#     if (group == "Lr34" & count == 21) {
#       plotter(chr, count, ylim = c(0, 1), smooth = FALSE, ab = "raw", yaxt = "s", xaxt = "s")
#       title("AMOVA Phi Stat between known Lr34 and lr34 Types, Chromosome 7D", line = -1)
#     } else if (group == "Lr22a" & count == 6) {
#       plotter(chr, count, ylim = c(0, 1), smooth = FALSE, ab = "raw", yaxt = "s", xaxt = "s")
#       title("AMOVA Phi Stat between known Lr22a and lr22a Types, Chromosome 2D", line = -1)
#     } else if (group == "Lr21" & count == 3) {
#       plotter(chr, count, ylim = c(0, 1), smooth = FALSE, ab = "raw", yaxt = "s", xaxt = "s")
#       title("AMOVA Phi Stat between known Lr21 and lr21 Types, Chromosome 1D", line = -1)
#     } else if (group == "Lr10" & count == 1) {
#       plotter(chr, count, ylim = c(0, 1), smooth = FALSE, ab = "raw", yaxt = "s", xaxt = "s")
#       title("AMOVA Phi Stat between known Lr10 and lr10 Types, Chromsome 1A", line = -1)
#     } else if (group == "Lr1" & count == 15) {
#       plotter(chr, count, ylim = c(0, 1), smooth = FALSE, ab = "raw", yaxt = "s", xaxt = "s")
#       title("AMOVA Phi Stat between known Lr1 and lr1 Types, Chromsome 5D", line = -1)
#     }
#     count <<- count + 1
#   })
#   title(xlab = "Marker by Chromosome Position", outer = T, cex.lab = 1.5, line = 1.5)
#   title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5, line = -1.5)
# }
# par(xpd = NA)
# par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
# plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# legend("bottom", legend = c("Raw 0.1% Threshold"), pch = 19, col = c("red"), cex = 0.6, horiz = T)
# dev.off()

