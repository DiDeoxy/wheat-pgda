source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "coord_flip", "element_text", "facet_grid", "geom_bar",
  "ggplot", "ggsave","ggtitle", "labs", "scale_fill_manual", "theme", "ylim"
)
import::from(magrittr, "%>%")
import::from(readr, "read_csv")

amova_table <- read_csv(
  file.path("results", "categorizations_AMOVA_variance.csv")
)

amova_table$Hierarchy <- factor(
  amova_table$Hierarchy,
  levels = c(
    "BP", "Era", "MTG", "MTG(BP): MTG", "MTG(BP): BP", "MTG(Era): MTG",
    "MTG(Era): Era", "MTG(BP(Era)): MTG", "MTG(BP(Era)): BP",
    "MTG(BP(Era)): Era", "MTG(Era(BP)): MTG", "MTG(Era(BP)): Era",
    "MTG(Era(BP)): BP", "BP-Era", "MTG(BP-Era): MTG", "MTG(BP-Era): BP-Era",
    "Clusters", "Clusters(BP-Era): Clusters", "Clusters(BP-Era): BP-Era"
  )
)

amova_table$type <- factor(
  c(
    "BP", "Era", "MTG", "MTG(BP)", "MTG(BP)", "MTG(Era)",
    "MTG(Era)", "MTG(BP(Era))", "MTG(BP(Era))", "MTG(BP(Era))",
    "MTG(Era(BP))", "MTG(Era(BP))", "MTG(Era(BP))", "BP-Era", 
    "MTG(BP-Era)", "MTG(BP-Era)", "Clusters",
    "Clusters(BP-Era)", "Clusters(BP-Era)"
  ),
  levels = c(
    "BP", "Era", "MTG", "MTG(BP)", "MTG(Era)", "MTG(BP(Era))", "MTG(Era(BP))",
    "BP-Era", "MTG(BP-Era)", "Clusters", "Clusters(BP-Era)"
  )
)

amova_table$class <- factor(
  c(
    "BP", "Era", "MTG", "MTG", "BP", "MTG", "Era", "MTG", "BP", "Era", "MTG",
    "Era", "BP", "BP-Era", "MTG", "BP-Era", "Clusters", "Clusters", "BP-Era"
  ),
  levels = c(
    "BP", "Era", "MTG", "BP-Era", "Clusters"
  )
)

max_var <- by(amova_table, amova_table$type, function (type_data) {
  sum(type_data$`% Variation`)
}) %>% max() %>% ceiling()

plots <- by(amova_table, amova_table$type, function (type_data) {
  if (type_data$type[1] == "MTG(BP(Era))") {
    type_data$class <- factor(
      type_data$class, levels = levels(type_data$class)[c(2, 1, 3:5)]
    )
    colour_set <- colour_set[c(2, 1, 3:5)]
  } else if (type_data$type[1] == "MTG(BP-Era)") {
    type_data$class <- factor(
      type_data$class, levels = levels(type_data$class)[c(4, 1:3, 5)]
    )
    colour_set <- colour_set[c(4, 1:3, 5)]
  }

  type_data %>%
    ggplot(aes(type, `% Variation`, fill = class)) +
      geom_bar(stat = "identity") +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 20)
      )  +
      ylim(0, max_var) +
      labs(fill = "Categorization") +
      scale_fill_manual(drop = FALSE, values = colour_set)
})

plots_matrix <- ggmatrix(
  plots, nrow = 1, ncol = 11,
  xlab = "Hierarchy",
  ylab = "% Variation",
  # title = "Genetic Distance Variance Partitioning",
  legend = c(1, 1)
) + theme(
  legend.position = "bottom", legend.box = "vertical",
  text = element_text(size = 20)
)

png(
  file.path("results", "categorization_AMOVA_variance.png"),
  width = 210, height = 180, units = "mm", res = 192
)
plots_matrix
dev.off()
