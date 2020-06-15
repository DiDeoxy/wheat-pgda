library(SNPRelate)
library(extrafont)

wheat <- snpgdsOpen("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Data\\Formatted\\wheat_phys_subset_sample.gds")
snp.pos <- read.gdsn(index.gdsn(wheat, "snp.position"))
snp.chrom <- factor(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
snp.id <- factor(read.gdsn(index.gdsn(wheat, "snp.id")))


png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\pca\\wheat_PCA_test1.png",
    family="Times New Roman", height = 11, width = 7, pointsize = 10, units = "in", res = 150)
par(mfrow = c(11,7), mar = c(0, 1, 1, 0))
for (ld in seq(0, 1, 0.1)) {
  for (bp in c(1e5, 1e6, 5e6, 1e7, 1e8, 1e9, 1e10)) {
    set.seed(1000)
    kept.snps.list <- snpgdsLDpruning(wheat, autosome.only = F, missing.rate = 0.1,
                                       ld.threshold = ld, slide.max.bp = bp)
    kept.snps <- unlist(kept.snps.list)
    pca <- snpgdsPCA(wheat, snp.id = kept.snps, num.thread = 4, autosome.only = F)
    pc.percent <- pca$varprop*100
    plot(pca$eigenvect[,1:2],
         pch = 19, cex = 0.8, xaxt = "n", yaxt = "n",
         main = paste("LD: ", ld, " BP: ", bp), cex.main = 0.75)
    mtext(side = 3, paste("SNPs: ", length(pca$snp.id)), line = -1, cex = 0.7, adj = 0)
    mtext(side = 2, paste(round(pc.percent[1], 2), "% ", round(pc.percent[2], 2), "%"), cex = 0.5)
  }
}
dev.off()

# png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\pca\\HRS_PCA_test2.png",
#     family="Times New Roman", width = 10, height = 10, pointsize = 10, units = "in", res = 150)
# par(mar = c(0, 1, 1, 0))
# pca <- snpgdsPCA(wheat, num.thread = 4, autosome.only = F)
# pc.percent <- pca$varprop*100
# plot(pca$eigenvect[,1:2],
#      pch = 19, cex = 0.8, xaxt = "n", yaxt = "n")
# mtext(side = 3, paste("SNPs: ", length(pca$snp.id)), line = -1, adj = 0)
# mtext(side = 2, paste(round(pc.percent[1], 2), "% ", round(pc.percent[2], 2), "%"))
# dev.off()

# png("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\Analysis\\Figures\\pca\\HRS_PCA_test3.png",
#     family="Times New Roman", width = 4, height = 6, pointsize = 10, units = "in", res = 150)
# par(mfrow = c(3,2), mar = c(0, 1, 1, 0))
# for (bp in c(1e6, 5e6, 1e7, 1e8, 1e9, 1e10)) {
#   set.seed(1000)
#   kept.snps.list <- snpgdsLDpruning(wheat, autosome.only = F, missing.rate = 0.1,
#                                     ld.threshold = 1, slide.max.bp = bp)
#   kept.snps <- unlist(kept.snps.list)
#   pca <- snpgdsPCA(wheat, snp.id = kept.snps, num.thread = 4, autosome.only = F)
#   pc.percent <- pca$varprop*100
#   plot(pca$eigenvect[,1:2],
#        pch = 19, cex = 0.8, xaxt = "n", yaxt = "n",
#        main = paste("LD: ", 1, " BP: ", bp), cex.main = 0.75)
#   mtext(side = 3, paste("SNPs: ", length(pca$snp.id)), line = -1, cex = 0.7, adj = 0)
#   mtext(side = 2, paste(round(pc.percent[1], 2), "% ", round(pc.percent[2], 2), "%"), cex = 0.5)
# }
# dev.off()
# snpgdsClose(wheat)

