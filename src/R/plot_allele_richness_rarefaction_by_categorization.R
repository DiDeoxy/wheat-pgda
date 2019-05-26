source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(
  ggplot2, "aes", "annotation_logticks", "element_blank", "element_text",
  "geom_line", "ggplot", "ggsave", "ggtitle","guide_legend", "guides", "labs",
  "scale_colour_manual", "scale_x_log10", "theme", "xlab", "ylab"
)
import::from(magrittr, "%>%")
import::from(parallel, "detectCores")
import::from(pgda, "allele_richness", "snpgds_parse")
import::from(readr, "read_rds")
import::from(scrime, "knncatimpute", "rowTables")
import::from(stringr, "str_c", "str_replace", "str_wrap")
import::from(tibble, "tibble")

wheat_data <- snpgds_parse(ld_gds)

geno_imputed <- wheat_data$geno %>%
    replace(. == 0, 1) %>%
    replace(. == 3, NA) %>%
    t() %>%
    knncatimpute() %>%
    t()

clusters <- factor(read_rds(hdbscan)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"
)

categorizations <- list(
  clusters, wheat_data$sample$annot$bp, wheat_data$sample$annot$era,
  wheat_data$sample$annot$mc, wheat_data$sample$annot$pheno,
  wheat_data$sample$annot$habit, wheat_data$sample$annot$colour,
  wheat_data$sample$annot$texture
)

categorization_names <- c(
  "HDBSCAN Clusters", "Breeding Program", "Era", "Market Class", "Phenotype",
  "Growth Habit", "Colour", "Texture"
)

categorization_colours <- list(
  colours_hdbscan_ar_legend, colours_bp, colours_era, colours_mc,
  colours_pheno_ar, colour_set[c(1, 22, 4)], colour_set[c(1, 22, 4)],
  colour_set[c(1, 4, 22)]
)

lapply(1:length(categorizations), function (i) {
  categorization <- categorizations[[i]]
  print(categorization_names[i])
  categories <- categorization %>% as.factor() %>% levels()

  rarefied <- lapply(seq_along(categories), function (j) {
    category <- categories[j]
    indivs <- which(categorization == category)
    if (length(indivs) > 1) {
      tibble(
        category = category, sample_size = 1:length(indivs),
        allele_richness = allele_richness(
          geno_imputed[, indivs], num_cores = detectCores()
        ) %>% colMeans()
      )
    }
  }) %>% do.call(rbind, .)

  plot <- rarefied %>% ggplot() +
    geom_line(
      aes(sample_size, allele_richness, colour = str_wrap(category, 10))
    ) +
    scale_x_log10() +
    scale_colour_manual(values = categorization_colours[[i]]) +
    ggtitle(
      str_c(
        "Rarefaction Curves of Average Allele\nRichness of Sub-samples By\n",
        categorization_names[i]
      )
    ) +
    xlab("Sample Size") +
    ylab("Average Allele Richness") +
    labs(colour = str_wrap(categorization_names[i], 5)) +
    guides(
      colour = guide_legend(nrow = (
        (rarefied$category %>% unique() %>% length()) / 3
      ) %>% ceiling())
    ) +
    theme(
      legend.position = "bottom", 
      legend.text = element_text(
        size = ifelse(length(categories) <= 4, 10,
          ifelse(length(categories) <= 8, 8, 4)
        )
      ),
      panel.grid.minor = element_blank()
    ) +
    annotation_logticks(sides = "b")

  ggsave(
    str_c(
      "allele_richness_rarefaction_by_",
      str_replace(categorization_names[i], " ", "_"), ".png"
    ),
    plot = plot, path = allele_richness, width = 100, height = 120, units = "mm"
  )
})
