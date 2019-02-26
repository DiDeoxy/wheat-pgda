# library(plyr)
# load the data from the gds object
wheat_data <- parse_gds(phys_gds)

# create a list of the ld between markers on each chromosome add NA column 
# and rows every 7.5 Mb so that gaps appear on image
wheat_gds <- snpgdsOpen(phys_gds)
ld_mat <- by(wheat_data$snp, wheat_data$snp$chrom, function (chrom) {
  gap_pos <- vector()
  prev <- 0
  dist <- 7.5
  for (cur in chrom$pos_mb) {
    if (cur - prev > dist) {
      gap_pos <- c(
        gap_pos,
        # finding the average positions between two markers that are
        # ~ dist apart
        cbind(seq(prev, cur, dist), rev(seq(cur, prev, -dist))) %>% rowMeans()
      )
    }
    prev <- cur
  }
  col_gaps_mat <- snpgdsLDMat(
      wheat_gds, method = "composite", snp.id = chrom$id, slide = -1
    )$LD %>%
    abs() %>%
    cbind(pos_mb = chrom$pos_mb) %>% 
    as.tibble() %>%
    full_join(tibble(pos_mb = gap_pos)) %>%
    arrange(pos_mb) %>%
    select(-pos_mb)
  both_gaps_mat <- t(col_gaps_mat) %>%
    cbind(pos_mb = chrom$pos_mb) %>%
    as.tibble() %>%
    full_join(tibble(pos_mb = gap_pos)) %>%
    arrange(pos_mb)
  list(pos_mb = both_gaps_mat$pos_mb,
       mat = both_gaps_mat %>% select(-pos_mb) %>% as.matrix())
})
snpgdsClose(wheat_gds)

# plot the matrices
cex <- 0.8
labels <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()
colours <- colorRampPalette(rev(brewer.pal(5, "RdYlBu")))(100)
png(
  file.path("results", "chrom_ld_heatmaps.png"), family = "Times New Roman",
  width = 135, height = 267, pointsize = 12, nits = "mm", res = 300
)
par(mfrow = c(7, 3), oma = c(7, 1, 3, 2), mar = c(0, 2, 0, 0))
for (chr in 1:21) {
  lim <- c(0, wheat_data$max_lengths[ifelse(chr %% 3, chr %% 3, 3)])
  if (chr %in% c(19, 20, 21)) {
    image(ld_mat[[chr]]$pos_mb, ld_mat[[chr]]$pos_mb, ld_mat[[chr]]$mat,
          col = colours, yaxt = "n", xaxt = "s", xlim = lim, ylim = lim)
    mtext(text = labels[chr], 2, line = 0.5, cex = cex)
  } else {
    image(ld_mat[[chr]]$pos_mb, ld_mat[[chr]]$pos_mb, ld_mat[[chr]]$mat,
          col = colours, yaxt = "n", xaxt = "n", xlim = lim, ylim = lim)
    mtext(text = labels[chr], 2, line = 0.5, cex = cex)
  }
}
title(xlab = "Marker Position in Mb", outer = T, cex.lab = 1.5, line = 2.5)
title(main = paste("Chromsomale LD Heatmaps"),
  outer = T, cex.main = 1.5, line = 1)
par(xpd = NA)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(-.70, -1.01,
  legend = c("Random Assortment", "Middling LD", "Complete LD"),
  fill = colours[c(1, 50, 100)], horiz = T)
dev.off()