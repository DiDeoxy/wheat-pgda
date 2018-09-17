# library(plyr)
library(tidyverse)
library(RColorBrewer)
library(extrafont)
library(SNPRelate)

# load data
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")

# create a list of the ld between markers on each chromosome add NA column 
# and rows every million bps so that gaps appear on image
full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_subset_sample.gds")
ld_mat <- by(data, snp_chrom, function (chrom) {
  gap_pos <- vector()
  prev = 0
  dist = 7.5
  for (cur in chrom$pos) {
    if (cur - prev > dist) {
      gap_pos <- c(gap_pos, 
                   # finding the average positions between two markers that are
                   # ~ dist apart
                   rowMeans(cbind(seq(prev, cur, dist), # forward
                                  rev(seq(cur, prev, -dist))))) # reverse
    }
    prev <- cur
  }
  col_gaps_mat <- snpgdsLDMat(full, method = "composite", snp.id = chrom$id,
    slide = -1)$LD %>%
    abs() %>%
    cbind(pos = chrom$pos) %>% 
    as.tibble() %>%
    full_join(tibble(pos = gap_pos)) %>%
    arrange(pos) %>%
    select(-pos) 
  both_gaps_mat <- t(col_gaps_mat) %>%
    cbind(pos = chrom$pos) %>%
    as.tibble() %>%
    full_join(tibble(pos = gap_pos)) %>%
    arrange(pos)
  list(pos = both_gaps_mat$pos,
       mat = both_gaps_mat %>% select(-pos) %>% as.matrix())
})
snpgdsClose(wheat)

# plot the matrices
cex <- 0.8
labels <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste, 
                            sep="")))
colours <- colorRampPalette(rev(brewer.pal(5, "RdYlBu")))(100)
png("Results\\loci\\LD\\full_LD_heatmap_7_5.png",
    family = "Times New Roman", width = 140, height = 287, pointsize = 12,
    units = "mm", res = 300)
par(mfrow = c(7,3), oma = c(7,1,3,2), mar = c(0, 2, 0, 0))
for (chr in 1:21) {
  lim = c(0, max_genome_lengths[ifelse(chr %% 3, chr %% 3, 3)])
  if (chr %in% c(19, 20, 21)) {
    image(ld_mat[[chr]]$pos, ld_mat[[chr]]$pos, ld_mat[[chr]]$mat,
          col = colours, yaxt = "n", xaxt = "s", xlim = lim, ylim = lim)
    mtext(text = labels[chr], 2, line = 0.5, cex = cex)
  } else {
    image(ld_mat[[chr]]$pos, ld_mat[[chr]]$pos, ld_mat[[chr]]$mat,
          col = colours, yaxt = "n", xaxt = "n", xlim = lim, ylim = lim)
    mtext(text = labels[chr], 2, line = 0.5, cex = cex)
  }
}
title(xlab = "Marker Position in Mb", outer = T, cex.lab = 1.5, line = 2.5)
title(main = paste("LD Heatmaps of Chromosmal Arms of All Varieties"),
      outer = T, cex.main = 1.5, line = 1)
par(xpd = NA)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(-.70, -1.01, 
       legend = c("Random Assortment", "Middling LD", "Complete LD"),
       fill = colours[c(1, 50, 100)], horiz = T)
dev.off()