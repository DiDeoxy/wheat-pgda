source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "full_join", "arrange", "select")
import::from(pgda, "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(magrittr, "%>%")
import::from(RColorBrewer, "brewer.pal")
import::from(SNPRelate, "snpgdsClose", "snpgdsLDMat", "snpgdsOpen")
import::from(tibble, "as_tibble", "tibble")

# library(plyr)
# load the data from the gds object
# wheat_data <- snpgds_parse(phys_gds)
wheat_data <- snpgds_parse(phys_gds)

max_lengths <- span_by_chrom(
  wheat_data$snp$chrom, wheat_data$snp$pos
) %>% max_lengths() / 1e6

# create a list of the ld between markers on each chromosome add NA column
# and rows every 7.5 Mb so that gaps appear on image
# wheat_gds <- snpgdsOpen(phys_gds)
wheat_gds <- snpgdsOpen(phys_gds)
ld_mat <- by(wheat_data$snp, wheat_data$snp$chrom, function (chrom) {
  gap_pos <- vector()
  prev <- 0
  dist <- 7.5e6
  for (cur in chrom$pos) {
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
    cbind(pos = chrom$pos) %>% 
    as_tibble() %>%
    full_join(tibble(pos = gap_pos)) %>%
    arrange(pos) %>%
    select(-pos)
  both_gaps_mat <- t(col_gaps_mat) %>%
    cbind(pos = chrom$pos) %>%
    as_tibble() %>%
    full_join(tibble(pos = gap_pos)) %>%
    arrange(pos)
  list(pos = both_gaps_mat$pos / 1e6,
       mat = both_gaps_mat %>% select(-pos) %>% as.matrix())
})
snpgdsClose(wheat_gds)

# plot the matrices
cex <- 0.8
labels <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()
# create a function for making a gradient of colours
colours <- colorRampPalette(colour_set[c(4, 2, 3, 5, 1)])(101)
png(
  file.path("results", "ld_heatmaps.png"), family = "Times New Roman",
  width = 120, height = 240, pointsize = 12, units = "mm", res = 300
)
par(mfrow = c(7, 3), oma = c(7, 1, 3, 2), mar = c(0, 2, 0, 0))
for (chr in 1:length(ld_mat)) {
  lim <- c(0, max_lengths[ifelse(chr %% 3, chr %% 3, 3)])
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
title(main = paste("LD Heatmaps by Chromosome"),
  outer = T, cex.main = 1.5, line = 1)
par(xpd = NA)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(-.60, -.98,
  legend = c("0", "0.25", "0.50", "0.75", "1"), title = "Abs. Composite LD",
  fill = colours[c(1, 26, 51, 76, 101)], horiz = T)
dev.off()

