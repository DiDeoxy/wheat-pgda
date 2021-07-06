# load file paths, colours, and functions
source("wheat-pgda/src/R/file_paths.R")
source("wheat-pgda/src/R/colour_sets.R")
import::from(gridExtra, "grid.arrange")
import::from(magrittr, "%>%")
import::from(pgda, "snpgds_parse")
import::from(
  ggplot2, "aes", "element_text", "ggplot", "geom_bar", "theme", "ylab", "xlab"
)
import::from(readr, "read_rds", "write_csv")
import::from(tibble, "as_tibble", "tibble")
import::from(tidyr, "gather")

# load the data from the gds object
all_data <- snpgds_parse(file.path(gds, "phys.gds"))
pruned_data <- snpgds_parse(ld_phys_gds)

class_data <- tibble(
  `Breeding Program` = all_data$sample$annot$bp %>% as.character(),
  Era = all_data$sample$annot$era %>% as.character(),
  `Major Trait Group` = all_data$sample$annot$mtg %>% as.character(),
  `Market Class` = all_data$sample$annot$mc %>% as.character()
) %>% gather(class, group)

plots <- lapply(class_data %>% split(class_data$class),
  function (class) {
    class %>%
      ggplot() +
      geom_bar(aes(group)) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
      ) +
      xlab(class[[1]]) +
      ylab("Count")
  }
)

# plot the matrix
png(
  file.path("results", "categoization_category_counts.png"),
  family = "Times New Roman", width = 1240, height = 1240
)
grid.arrange(
  grobs = plots, nrow = 2, ncol = 2
)
dev.off()

genos_names <- t(all_data$genotypes) %>% replace(3, NA)
colnames(genos_names) <- all_data$snp$id

clusters <- factor(read_rds(hdbscan_rds)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"
)
clusters <- clusters %>% as.character()
removed <- which(! all_data$sample$id %in% pruned_data$sample$id)
indices <- c(1:length(clusters), removed + 0.5)
clusters <- c(clusters, rep("N/A", length(removed)))[order(indices)]

write_csv(
  sapply(all_data$sample$annot, function(annot) as.character(annot)) %>%
    cbind(
      cultivar = all_data$sample$id, ., cluster = clusters, genos_names
    ) %>% as_tibble(),
  file.path("results", "cultivars_and_metadata.csv")
)