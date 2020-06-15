library(SNPRelate)
library(RColorBrewer)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")
source("Analysis\\R\\functions\\colour_sets.R")

load("Data\\Intermediate\\CCP\\full_kmeans.RDATA")
groups <- c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6")

for (maf in c(seq(0, 0.2, 0.01), seq(.25, .45, .05))) {
  ## positive
  ## setting up the data
  gdsSubset <- paste0("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_", maf*100,".gds")
  source("Analysis\\R\\functions\\data_loading.R")

  ## PCA
  wheat <- snpgdsOpen(paste0("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_", maf*100,".gds"))
  pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
  snpgdsClose(wheat)
  pc.percent <- pca$varprop*100

  png(paste0("Results\\pca\\maf\\pca_maf_", maf*100,"_to_50.png"), family="Times New Roman",
      width = 6, height = 7, pointsize = 10, units = "in", res = 300)
  lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
  pairs(pca$eigenvect[,1:3],
        labels=lbls, pch = 19, cex = 0.8,
        main = paste("MAF", maf, "to 0.50"),
        col = coloursIcl[bests],
        oma = c(5,2,5,2)
  )
  par(fig = c(0, 1, 0, 1), oma = c(0.5, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = groups, pch=19, cex = 0.6, horiz = T, col = coloursDend)
  text(-0.7, 0.95, labels = paste("Number of markers:", length(snp.id)))
  dev.off()

  if (maf != 0.00) {
    ## minus
    gdsSubset <- paste0("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_minus_", maf*100,".gds")
    source("Analysis\\R\\functions\\data_loading.R")

    ## designation popcode
    # pop.code <- as.factor(replace(desig, desig == "N/A", "UNKNOWN"))

    ## PCA
    wheat <- snpgdsOpen(paste0("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_minus_", maf*100,".gds"))
    pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
    snpgdsClose(wheat)
    pc.percent <- pca$varprop*100

    png(paste0("Results\\pca\\maf\\pca_maf_0_to_", maf*100,".png"), family="Times New Roman",
        width = 6, height = 6.5, pointsize = 10, units = "in", res = 300)
    lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
    pairs(pca$eigenvect[,1:3],
          labels=lbls, pch = 19, cex = 0.8,
          main = paste("MAF 0 to", maf),
          col = coloursIcl[bests],
          oma = c(5,2,5,2)
    )
    par(fig = c(0, 1, 0, 1), oma = c(0.5, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = groups,
           pch=19, cex = 0.6, horiz = T, col = coloursDend)
    text(-0.7, 0.95, labels = paste("Number of markers:", length(snp.id)))
    dev.off()
  }
  
  gdsSubset <- paste0("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_slice_", maf*100, "_to_", (maf+0.05)*100, ".gds")
  source("Analysis\\R\\functions\\data_loading.R")
 
  ## PCA
  wheat <- snpgdsOpen(paste0("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_slice_", maf*100, "_to_", (maf+0.05)*100, ".gds"))
  pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
  snpgdsClose(wheat)
  pc.percent <- pca$varprop*100
  
  png(paste0("Results\\pca\\maf\\pca_maf_slice_", maf*100, "_to_", (maf+0.05)*100, ".png"), family="Times New Roman",
      width = 6, height = 7, pointsize = 10, units = "in", res = 300)
  lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
  pairs(pca$eigenvect[,1:3],
        labels=lbls, pch = 19, cex = 0.8,
        main = paste("MAF Slice", maf, "to", (maf+0.05)),
        col = coloursIcl[bests],
        # col = coloursDesig[pop.code],
        oma = c(5,2,5,2)
  )
  par(xpd = T)
  par(fig = c(0, 1, 0, 1), oma = c(0.5, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = groups,
         pch=19, cex = 0.6, horiz = T, col = coloursDend)
  text(-0.7, 0.95, labels = paste("Number of markers:", length(snp.id)))
  dev.off()
}