source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "bind_rows")
import::from(emdbook, "lseq")
import::from(GeneticSubsetter, "PicCalc")
import::from(
  ggplot2, "aes", "element_text", "geom_point", "geom_smooth", "ggplot",
  "ggsave", "ggtitle","guide_legend", "guides", "labs", "scale_colour_manual",
  "theme", "xlab", "ylab"
)
import::from(hierfstat, "allelic.richness")
import::from(magrittr, "%>%")
import::from(parallel, "detectCores", "mclapply")
import::from(pgda, "calc_eh", "snpgds_parse")
import::from(readr, "read_rds")
import::from(stringr, "str_c")
import::from(tibble, "add_column", "add_row", "as_tibble", "tibble")

wheat_data <- snpgds_parse(phys_gds)

wheat_ar <- wheat_data$genotypes %>%
  replace(. == 0, 1) %>%
  replace(. == 3, NA) %>%
  t() %>%
  as.data.frame()

# wheat_test <- wheat_data$genotypes %>%
#   replace(. == 0, 1) %>%
#   replace(. == 3, NA) %>%
#   t() %>%
#   as.data.frame() %>%
#   add_column(clusters %>% as.character(), .before = 1)
# test <- allelic.richness(wheat_test)
# test$Ar %>% head()
# test$Ar %>% colMeans() %>% str()

clusters <- factor(read_rds(hdbscan)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1 (HRW)", "Cluster 2 (SWS)", "Cluster 3 (CWES)",
  "Cluster 4 (CPSR/W)", "Cluster 5 (HRS)"
)

categorizations <- list(
  clusters, wheat_data$sample$annot$bp, wheat_data$sample$annot$era,
  wheat_data$sample$annot$mc, wheat_data$sample$annot$pheno,
  wheat_data$sample$annot$habit, wheat_data$sample$annot$colour,
  wheat_data$sample$annot$texture
)
categorization_names <- c(
  "clusters", "breeding_programs", "era", "market_class", "phenotype",
  "growth_habit", "colour", "texture"
)
categorization_colours <- list(
  colours_hdbscan_legend, colours_bp, colours_era, colours_mc, colours_pheno,
  colour_set[c(1, 22, 4)], colour_set[c(1, 22, 4)], colour_set[c(1, 22, 4)]
)

lapply(1, function (i) {
  subset_sizes <- lseq(2, 60, 15) %>% floor() %>% unique()
  categorization <- categorizations[[i]]
  category_sizes <- categorization %>% table()
  categories <- categorization %>% as.factor() %>% levels()

  rarefied <- mclapply(subset_sizes, function (subset_size) {
    cats_of_size <- categories[which(category_sizes >= subset_size)]
    kept_indivs <- which(categorization %in% cats_of_size)
    temp <- wheat_ar[kept_indivs, ] %>%
      add_column(
        categorization = categorization[kept_indivs] %>% as.character(),
        .before = 1
      )

    print(str_c(categorization_names[i], ": ", subset_size))
    tibble(
      category = temp$categorization %>% as.factor() %>% levels(),
      subset_size = subset_size,
      pic = allelic.richness(temp, min.n = subset_size)$Ar %>% colMeans()
    )
  }, mc.cores = detectCores()) %>% do.call(rbind, .)

  plot <- rarefied %>% ggplot() +
    geom_point(
      aes(subset_size, pic, colour = category), size = 0.5, alpha = 0.2,
      pch = 16
    ) +
    geom_smooth(aes(jitter(subset_size), pic, colour = category), se = FALSE) +
    scale_colour_manual(values = categorization_colours[[i]]) +
    labs(colour = "Category") +
    ggtitle(
      str_c(
        "Rarefaction Curves of Average\nAllelic Richness of Sub-samples\nBy ",
        categorization_names[i]
      )
    ) +
    xlab("Sample Size") +
    ylab("Average Allelic Richness") +
    guides(colour = guide_legend(nrow = 3)) +
    theme(
      legend.position = "bottom", legend.text = element_text(size = 4)
    )
  ggsave(
    str_c("Allelic_richness_rarefaction_by_", categorization_names[i], ".png"), plot = plot, 
    path = PIC, width = 100, height = 120, units = "mm"
  )
})