# load file paths, colours, and functions
source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(gridExtra, "grid.arrange")
import::from(magrittr, "%>%")
import::from(pgda, "snpgds_parse")
import::from(
  ggplot2, "aes", "element_text", "ggplot", "geom_bar", "theme", "ylab", "xlab"
)
import::from(readr, "write_csv")
import::from(tibble, "tibble")
import::from(tidyr, "gather")

# load the data from the gds object
wheat_data <- snpgds_parse(file.path(gds, "phys.gds"))

class_data <- tibble(
  `Breeding Program` = wheat_data$sample$annot$bp %>% as.character(),
  Era = wheat_data$sample$annot$era %>% as.character(),
  `Major Trait Group` = wheat_data$sample$annot$mtg %>% as.character(),
  `Market Class` = wheat_data$sample$annot$mc %>% as.character()
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
  family = "Times New Roman", width = 160, height = 192, pointsize = 5,
  units = "mm", res = 192
)
grid.arrange(
  grobs = plots, nrow = 2, ncol = 2
)
dev.off()

genos_names <- t(wheat_data$genotypes) %>% replace(3, NA)
colnames(genos_names) <- wheat_data$snp$id

write_csv(
  sapply(wheat_data$sample$annot, function(annot) as.character(annot)) %>%
    cbind(cultivar = wheat_data$sample$id, ., genos_names) %>%
    as.data.frame(),
  file.path("results", "cultivars_and_metadata.csv")
)