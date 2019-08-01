# load file paths, colours, and functions
source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(gridExtra, "grid.arrange")
import::from(magrittr, "%>%")
import::from(pgda, "snpgds_parse")
import::from(
  ggplot2, "aes", "element_text", "ggplot", "geom_bar", "theme", "ylab", "xlab"
)
import::from(tibble, "as_tibble")

# load the data from the gds object
wheat_data <- snpgds_parse(file.path(gds, "phys.gds"))

class_data <- rbind (
  cbind("Breeding Program", wheat_data$sample$annot$bp %>% as.character()),
  cbind("Era", wheat_data$sample$annot$era %>% as.character()),
  cbind("Major Trait Group", wheat_data$sample$annot$pheno %>% as.character()),
  cbind("Market Class", wheat_data$sample$annot$mc %>% as.character())
) %>% as_tibble()
rownames(class_data) <- c()
colnames(class_data) <- c("Class", "Group")

table(class_data)

plots <- lapply(class_data %>% split(class_data$Class),
  function (class) {
    class %>%
      ggplot() +
      geom_bar(aes(Group)) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)
      ) +
      xlab(class$Class[1]) +
      ylab("Count")
  }
)

# plot the matrix
png(
  file.path("results", "test.png"),
  family = "Times New Roman", width = 160, height = 192, pointsize = 5,
  units = "mm", res = 192
)
grid.arrange(
  grobs = plots, nrow = 2, ncol = 2
)
dev.off()
