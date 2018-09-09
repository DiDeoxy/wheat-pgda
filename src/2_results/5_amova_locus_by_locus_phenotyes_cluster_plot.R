library(extrafont)
library(SNPRelate)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA")
source("src\\R_functions\\funcs_draw_loci.R")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")
snp.pos <- snp.pos/1000000

# aligned genes
mainGenes <- get(load("Data\\Intermediate\\Aligned_genes\\top_main_genes_contigs.RData"))
resiGenes <- get(load("Data\\Intermediate\\Aligned_genes\\top_resistance_genes_contigs.RData"))
knownGenes <- get(load("Data\\Intermediate\\Aligned_genes\\known_genes_groups.RData"))
mainGenes$Pos <- mainGenes$Pos/1000000
resiGenes$Pos <- resiGenes$Pos/1000000
knownGenes$Pos <- knownGenes$Pos/1000000

main <- by(mainGenes, snp.chrom, function (x) x)
head(main[[6]])
resi <- by(resiGenes, snp.chrom, function (x) x)
head(resi[[6]])

# stuff <- by(mainGenes, snp.chrom, function (x) {x})
# stuff[[11]]

gene.pos <- 16.3566
for (group in c("HRS_HRW", "HRS_SWS", "HRW_SWS")[1]) {
  
  amova <- get(load(paste0("Data\\Intermediate\\Adegenet\\amova_", group, ".RData")))
  phiAll <- phi.all(amova)
  signif <- quantile(phiAll, prob = 0.975, na.rm = T)
  snp.pos[which(phiAll < signif)] <- NA
  
  blah <- by(data.frame(snp.pos, snp.id), snp.chrom, function (chr) {
    chr[,1]
    # print(chr[,1] - gene.pos)
    # print(min(abs(chr[,1] - gene.pos)))
    # return(as.character(chr[which.min(abs(chr[,1] - gene.pos)),2]))
  })
  print(blah)
}

snp.pos[which(snp.id == "wsnp_Ex_rep_c66689_65010988")] - gene.pos

snps <- list()
for (group in c("HRS_SWS", "HRS_HRW", "HRW_SWS")) {
  
  amova <- get(load(paste0("Data\\Intermediate\\Adegenet\\amova_", group, ".RData")))
  phiAll <- phi.all(amova)
  signif <- quantile(phiAll, prob = 0.975, na.rm = T)
  
  # snps[[group]] <- by(data.frame(snp.id, phiAll, snp.pos), snp.chrom, function (chr) {
  #   return(chr[which(chr[,2] > signif),3])
  # })

  # phiStats <- by(phiAll, snp.chrom, function (chr) {
  #   return(cbind(mean(abs(diff(chr))), sd(abs(diff(chr)))))
  # })
  # order <- c(seq(1, 19, 3), seq(2, 20, 3), seq(3, 21, 3), 22)
  # labels <- c(as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep=""))), "All")
  # phiTable <- data.frame()
  # phiTable <- do.call("rbind", phiStats)
  # phiTable <- rbind(phiTable, cbind(mean(phiAll), sd(phiAll)))
  # row.names(phiTable) <- labels
  # colnames(phiTable) <- c("Mean", "Standard Deviation")
  # phiTable <- phiTable[order,]
  # write.csv(phiTable, file = paste0("Data\\Intermediate\\phi stats\\", group, ".csv"))
  
  ## plotting
  png(paste0("Results\\loci\\amova\\full_amova_", group, ".png"),
      family = "Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
  plots(snp.pos, phiAll, mainGenes, resiGenes, snp.chrom, ylim = c(signif, 1), oma = c(4,0,3,4), plotGenes = T)
  title(xlab = "Marker Position in Mb", outer = T, cex.lab = 1.5, line = 2.5)
  mtext("AMOVA Phi Stat", outer = T, cex.lab = 1.5, line = 2.5, side = 4)
  names <- unlist(strsplit(group, "_"))
  if (names[1] == "HRS") names[1] <- "CHRS"

  par(xpd = NA)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  title(paste("Top 2.5% AMOVA Phi Statistics By Locus Between", names[1], "and", names[2], "Wheat Varieties"),
        outer = T, cex.main = 1.5, line = -1.5)
  legend(-1.01, -1.01, legend = c("Phenotype and Dwarfing Genes", "Resistance Genes"),
         pch = 19, col = c(colourSet[15], colourSet[19]), cex = 0.6, horiz = T)
  dev.off()
}
snps

# 
png(paste0("Results\\loci\\amova\\full_amova_genes.png"), family = "Times New Roman", width = 6, height = 6, pointsize = 10, units = "in", res = 300)
par(mfrow = c(5,1), mar = c(1,4,1,2), oma = c(4,0,0,0))
emptyFrame <- data.frame(one = rep(NA, 15572), two = rep(NA, 15572), three = rep(NA, 15572))
for (group in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
  amova <- get(load(paste0("Data\\Intermediate\\Adegenet\\amova_", group, ".RData")))
  phiAll <- phi.all(amova)
  signif <- quantile(phiAll, prob = 0.975, na.rm = T)

  ## plotting
  count <- 1
  by(cbind(snp.pos, phiAll, emptyFrame, knownGenes), snp.chrom, function (chr) {
    if (group == "Lr34" & count == 21) {
      print("Lr34")
      print(chr[,7])
      print(chr[,1][which(chr[,2] >= signif)])
      print(chr[,2][which(chr[,2] >= signif)])
      plotter(chr, count, ylim = c(signif, 1), genes = TRUE, yaxt = "s", xaxt = "s", colour = TRUE)
      title("Top 2.5% of AMOVA Phi Stat between known Lr34 and lr34 Types, Chromosome 7D", line = -1)
    } else if (group == "Lr22a" & count == 6) {
      print("Lr22a")
      print(chr[,7])
      print(chr[,1][which(chr[,2] >= signif)])
      print(chr[,2][which(chr[,2] >= signif)])
      plotter(chr, count, ylim = c(signif, 1), genes = TRUE, yaxt = "s", xaxt = "s", colour = TRUE)
      title("Top 2.5% of AMOVA Phi Stat between known Lr22a and lr22a Types, Chromosome 2D", line = -1)
    } else if (group == "Lr21" & count == 3) {
      print("Lr21")
      print(chr[,7])
      print(chr[,1][which(chr[,2] >= signif)])
      print(chr[,2][which(chr[,2] >= signif)])
      plotter(chr, count, ylim = c(signif, 1), genes = TRUE, yaxt = "s", xaxt = "s", colour = TRUE)
      title("Top 2.5% of AMOVA Phi Stat between known Lr21 and lr21 Types, Chromosome 1D", line = -1)
    } else if (group == "Lr10" & count == 1) {
      print("Lr10")
      print(chr[,7])
      print(chr[,1][which(chr[,2] >= signif)])
      print(chr[,2][which(chr[,2] >= signif)])
      plotter(chr, count, ylim = c(signif, 1), genes = TRUE, yaxt = "s", xaxt = "s", colour = TRUE)
      title("Top 2.5% of AMOVA Phi Stat between known Lr10 and lr10 Types, Chromsome 1A", line = -1)
    } else if (group == "Lr1" & count == 15) {
      print("Lr1")
      print(chr[,7])
      print(chr[,1][which(chr[,2] >= signif)])
      print(chr[,2][which(chr[,2] >= signif)])
      plotter(chr, count, ylim = c(signif, 1), genes = TRUE, yaxt = "s", xaxt = "s", colour = TRUE)
      title("Top 2.5% of AMOVA Phi Stat between known Lr1 and lr1 Types, Chromsome 5D", line = -1)
    }
    count <<- count + 1
  })
  title(xlab = "Marker by Chromosome Position in Mb", outer = T, cex.lab = 1.5, line = 1.5)
  title(ylab = "AMOVA Phi Statistic", outer = T, cex.lab = 1.5, line = -1.5)
}
par(xpd = NA)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
dev.off()

