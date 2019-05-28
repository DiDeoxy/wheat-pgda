# load file paths and functions
source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(
  circlize, "circos.initialize", "circos.clear", "circos.dendrogram", 
  "circos.par", "circos.text", "circos.track"
)
import::from(dendextend, "color_branches", "set")
import::from(magrittr, "%>%")
import::from(pgda, "draw_rects", "snpgds_parse")
import::from(readr, "read_rds")
import::from(SNPRelate, "snpgdsClose", "snpgdsIBS", "snpgdsOpen")
import::from(stringr, "str_c")

## setting up the data
wheat_data <- snpgds_parse(ld_gds)

clusters <- factor(read_rds(hdbscan)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"
)

table(wheat_data$sample$annot$mc)

wheat_gds <- snpgdsOpen(ld_gds)
## making the distance object
ibs_dist <- as.dist(1 - snpgdsIBS(wheat_gds, autosome.only = F)$ibs)
snpgdsClose(wheat_gds)

upgma_dend <- hclust(ibs_dist) %>%
  as.dendrogram(method = "average") %>%
  color_branches(k = 6, col = colours_dend) %>%
  set("branches_lwd", 1.5)
label_order <- order.dendrogram(upgma_dend)

## drawing the circos plot
png(
  file.path("results", "dend.png"), family = "Times New Roman", width = 210,
  height = 210, pointsize = 15, units = "mm", res = 500
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
  wheat_data$sample$annot$pheno, colours_pheno, label_order, colors()[525]
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
  "topleft", legend = levels(wheat_data$sample$annot$pheno), box.lwd = 2,
  title = "Phenotype", pch = pch, col = colours_pheno, cex = cex,
  bg = colors()[525]
)
legend(
  "center", legend = levels(clusters), title = "HDBSCAN Clusters", box.lwd = 2,
  pch = pch, col = colours_hdbscan_legend, cex = cex, bg = colors()[109]
)

title(
  main = str_c(
    "UPGMA Dendrogram of 365 Varieties\nwith HDBSCAN Clusters and ",
    "Categorical Data In Surrounding Rows"
  ),
  cex.main = 0.7
)
dev.off()

################################################################################
# out put some table summarizing the groupings
clusters <- read_rds(hdbscan)$cluster %>% replace(. == 0, "Noise") 

table(data.frame(wheat_data$sample$annot$pheno, clusters)) %>%
  as.data.frame.matrix() %>%
  write.csv(file.path("results", "pheno_cluster_table.csv"), quote = FALSE)
table(data.frame(wheat_data$sample$annot$mc, clusters)) %>%
  as.data.frame.matrix() %>%
  write.csv(file.path("results", "mc_cluster_table.csv"), quote = FALSE)