# library(plyr)
# library(tidyverse)
# library(reshape2)
# library(GGally)
# # install.packages("raster")
# library(raster)
# # install.packages("rasterVis")
# library(rasterVis)
# library(ggrepel)
# library(extrafont)
# library(SNPRelate)

# # load data
# gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
# source("src\\R_functions\\data_loading.R")

# # size of a megabase, used to divide the bp positions
# mb <- 1000000

# # find the max position of any marker on each genome for xlims
# chrom_lengths <- by(snp_pos / mb, snp_chrom, max)
# max_genome_lengths <- data.frame(A = max(chrom_lengths[seq(1, 19, 3)]),
#                                  B = max(chrom_lengths[seq(2, 20, 3)]),
#                                  D = max(chrom_lengths[seq(3, 21, 3)]))

# # create a list of the ld between markers on each chromosome
# full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_subset_sample.gds")
# plots <- list()
# ld_mat <- tibble()
# for (i in 1:21) {
#   pos <- snp_pos[which(snp_chrom == i)]
#   ld_mat <- abs(snpgdsLDMat(
#     full, method = "composite", snp.id = snp_id[which(snp_chrom == i)],
#     slide = -1)$LD) %>% cbind(pos) %>% as.tibble() %>%
#     full_join(tibble(pos = seq(0, max(pos), 1e6))) %>% 
#     dplyr::select(-pos)
#   ld_mat <- t(ld_mat) %>%
#     cbind(pos) %>%
#     as.tibble() %>%
#     full_join(tibble(pos = seq(0, max(pos), 1e6))) %>%
#     dplyr::select(-pos)
#   # plots[[i]] <- raster(ld_mat, xmn = 0, xmx = chrom_lengths[i], ymn = 0,
#   #                      ymx = chrom_lengths[i]) %>%
#   #                 gplot() +
#   #                   geom_raster(aes(fill = value))
#   break
# }
# snpgdsClose(wheat)
# ?raster
# as.vector(ld_mat[1,1200:1205])
# str(ld_mat)
# dim(ld_mat)
# ld_mat[, 1068]
# ?seq
# seq(5, 1e6, 1e7)
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


