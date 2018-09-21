library(tidyverse)
library(magrittr)
library(GGally)
library(SNPRelate)
library(extrafont)
# install.packages("rgl")
# library(rgl)

source("src\\R_functions\\colour_sets.R")

gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds"
source("src\\R_functions\\data_loading.R")

## PCA
full <- snpgdsOpen(
  "Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds"
)
pca <- snpgdsPCA(full, num.thread = 4, autosome.only = F)
snpgdsClose(full)

pc.percent <- pca$varprop * 100

cluster <- read_rds("Data\\Intermediate\\dbscan\\full_hdbscan.rds")$cluster

ggpairs_pca <- pca$eigenvect[, 1:3] %>%
  cbind(Cluster = cluster, .) %>%
  as.tibble() %>%
  rename(PC1 = V1, PC2 = V2, PC3 = V3) %<>%
  mutate(Cluster = factor(cluster)) %>%
  ggpairs(
    mapping = aes(colour = Cluster),
    title = paste(
      "First Three PCs With Varieties Coloured by",
      "HDBSCAN Clusters"
    )
  )
for (i in 1:ggpairs_pca$nrow) {
  for (j in 1:ggpairs_pca$ncol) {
    ggpairs_pca[i, j] <- ggpairs_pca[i, j] +
      scale_fill_manual(values = colours_dbscan) +
      scale_color_manual(values = colours_dbscan)
  }
}
ggpairs_pca[4, 1] <- ggpairs_pca[4, 1] +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous(breaks = c(0, 45))

png("Results\\pca\\full_pca_dbscan.png",
  family = "Times New Roman",
  width = 200, height = 200, pointsize = 5, units = "mm", res = 300
)
ggpairs_pca
dev.off()

# # 3D plot
# plot3d(pca$eigenvect[,1:3], col = colours_dbscan[fullHdbscan$cluster + 1],
#        type = "s", size = 1, xlab = "Eigenvector 1", ylab = "Eigenvector 2",
#        zlab = "Eigenvector 3")