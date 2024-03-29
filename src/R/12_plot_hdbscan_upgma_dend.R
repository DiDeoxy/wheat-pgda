# load file paths and functions
source("wheat-pgda/src/R/file_paths.R")
source("wheat-pgda/src/R/colour_sets.R")
import::from(
  circlize, "circos.initialize", "circos.clear", "circos.dendrogram", 
  "circos.par", "circos.text", "circos.track"
)
import::from(dendextend, "color_branches", "set")
import::from(magrittr, "%>%")
import::from(pgda, "draw_rects", "snpgds_parse")
import::from(readr, "read_rds", "write_csv")
import::from(SNPRelate, "snpgdsClose", "snpgdsIBS", "snpgdsOpen")
import::from(stringr, "str_c")

## setting up the data
wheat_data <- snpgds_parse(ld_phys_gds)

clusters <- factor(read_rds(hdbscan_rds)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"
)

# write_csv(
#   matrix(wheat_data$sample$id[which(clusters == "Noise")], 7, 5) %>%
#     as.data.frame(),
#   path = file.path("results", "noise_cultivars.csv")
# )

wheat_gds <- snpgdsOpen(ld_phys_gds)
## making the distance object
ibs_dist <- as.dist(1 - snpgdsIBS(wheat_gds, autosome.only = FALSE)$ibs)
snpgdsClose(wheat_gds)

upgma_dend <- hclust(ibs_dist) %>%
  as.dendrogram(method = "average") %>%
  color_branches(k = 9, col = colours_dend) %>%
  set("branches_lwd", 5)
label_order <- order.dendrogram(upgma_dend)

## drawing the circos plot
png(
  file.path("results", "dend.png"), family = "Times New Roman", width = 2480,
  height = 2480, pointsize = 60
)
circos.par(
  cell.padding = c(0, 0, 0, 0), gap.degree = 0.5, track.margin = c(0.005, 0.005)
)
circos.initialize("foo", xlim = c(0, length(label_order)), sector.width = 1)
circos.track(
  ylim = c(0, 1), track.height = 0.15, bg.border = NA,
  panel.fun = function(x, y) {
    circos.text(
      1:length(label_order), rep(0, length(label_order)),
      wheat_data$sample$id[label_order], facing = "clockwise",
      niceFacing = TRUE, cex = 0.25, adj = c(0, -0.2), font = 2
    )
  }
)

draw_rects(
  wheat_data$sample$annot$era, colours_era, label_order, colors()[468]
)
draw_rects(
  wheat_data$sample$annot$bp, colours_bp, label_order, colors()[383]
)
draw_rects(
  wheat_data$sample$annot$mc, colours_mc, label_order, colors()[105]
)
draw_rects(
  wheat_data$sample$annot$mtg, colours_mtg, label_order, colors()[525]
)
draw_rects(
  clusters, colours_hdbscan, label_order, colors()[109]
)

max_height <- max(attr(upgma_dend, "height"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
  circos.dendrogram(upgma_dend, max_height = max_height)
}, track.height = 0.3, bg.border = NA)

circos.clear()

## Legends
pch <- 19
cex <- 0.45
legend(
  "topright", legend = levels(wheat_data$sample$annot$era), box.lwd = 2,
  title = "Period of Release", pch = pch, col = colours_era, cex = cex,
  bg = colors()[468]
)
legend(
  "bottomright", legend = levels(wheat_data$sample$annot$bp), box.lwd = 2,
  title = "Breeding Program/Origin", pch = pch, col = colours_bp, cex = 0.35,
  bg = colors()[383]
)
legend(
  "bottomleft", legend = levels(wheat_data$sample$annot$mc), box.lwd = 2,
  title = "Market Class", pch = pch, col = colours_mc, cex = cex,
  bg = colors()[105]
)
legend(
  "topleft", legend = levels(wheat_data$sample$annot$mtg), box.lwd = 2,
  title = "Major Trait Group", pch = pch, col = colours_mtg, cex = cex,
  bg = colors()[525]
)
legend(
  "center", legend = levels(clusters), title = "HDBSCAN Clusters", box.lwd = 2,
  pch = pch, col = colours_hdbscan_legend, cex = cex, bg = colors()[109]
)

title(
  main = str_c(
    "Clustering by LD Pruned Markers"
    # "UPGMA Dendrogram of 365 Varieties\nwith HDBSCAN Clusters and ",
    # "Categorical Data In Surrounding Rows"
  ),
  cex.main = 0.7
)

dev.off()

################################################################################
# out put some table summarizing the groupings
clusters <- read_rds(hdbscan_rds)$cluster %>% replace(. == 0, "Noise") 

table(data.frame(wheat_data$sample$annot$mtg, clusters)) %>%
  as.data.frame.matrix() %>%
  write.csv(file.path("results", "mtg_cluster_table.csv"), quote = FALSE)
table(data.frame(wheat_data$sample$annot$mc, clusters)) %>%
  as.data.frame.matrix() %>%
  write.csv(file.path("results", "mc_cluster_table.csv"), quote = FALSE)