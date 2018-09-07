library(SNPRelate)
library(extrafont)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

gdsSubset <- paste0("Data\\Intermediate\\GDS\\wheat_phys_subset_sample.gds")
source("src\\functions\\data_loading.R")
snp.pos <- snp.pos/1000000

# making list of ld heatmaps
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, sep="")))
# "wheat", "HRS", "HRW", "SWS"
for (subset in c("wheat", "HRS", "HRW", "SWS")) {
  wheat <- snpgdsOpen(paste0("Data\\Intermediate\\GDS\\", subset, "_phys_subset_sample.gds"))
  ld_mat <- list()
  for(i in levels(as.factor(snp.chrom))) {
    ld_mat[[labels[as.numeric(i)]]] <- c(snpgdsLDMat(wheat, method = "composite",
                                                     snp.id = snp.id[which(snp.chrom == i)], 
                                                     slide = -1, verbose = F),
                                         list(pos = snp.pos[which(snp.chrom == i)]))
  }
  snpgdsClose(wheat)

  # plotting function
  cex <- 0.8
  pch <- 20
  xlim = c(0, max(snp.pos))
  ylim = c(0, max(snp.pos))
  colours <- colorRampPalette(c("khaki","green","green4","violet","purple"))(100)

  png(paste0("Results\\loci\\LD\\", subset, "_LD_heatmap.png"),
      family="Times New Roman", width = 3.5, height = 7, pointsize = 10, units = "in", res = 300)
  par(mfrow = c(7,3), oma = c(7,1,3,2), mar = c(0, 2, 0, 0))
  count <- 1
  for (name in names(ld_mat)) {
    if (count %in% c(19, 20, 21)) {
      image(ld_mat[[name]]$pos, ld_mat[[name]]$pos, ld_mat[[name]]$LD^2, col = colours,
            yaxt="n", xaxt = "s", xlim = xlim, ylim = ylim)
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    } else {
      image(ld_mat[[name]]$pos, ld_mat[[name]]$pos, ld_mat[[name]]$LD^2, col = colours,
            yaxt="n", xaxt = "n", xlim = xlim, ylim = ylim)
      mtext(text = labels[count], 2, line = 0.5, cex = cex)
    }
    count = count + 1
  }
  title(xlab = "Marker Position in Mb", outer = T, cex.lab = 1.5, line = 2.5)
  if (subset == "wheat") {
    title(main = paste("LD Heatmaps of Chromosmal Arms of All Varieties"), outer = T, cex.main = 1.5, line = 1)
  } else {
    title(main = paste("LD Heatmaps of Chromosmal Arms of", subset, "Varieties"), outer = T, cex.main = 1.5, line = 1)
  }

  par(xpd = NA)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(-.92, -.99, legend = c("Random Assortment", "Middling LD", "Complete LD"), fill = colours[c(1, 50, 100)], horiz = T)
  dev.off()
}

