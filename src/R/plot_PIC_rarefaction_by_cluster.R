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
import::from(magrittr, "%>%")
import::from(parallel, "detectCores", "mclapply")
import::from(pgda, "calc_eh", "snpgds_parse")
import::from(readr, "read_rds")
import::from(stringr, "str_c", "str_replace", "str_wrap")
import::from(tibble, "add_column", "add_row", "as_tibble", "tibble")

wheat_data <- snpgds_parse(phys_gds)

wheat_formatted <- wheat_data$genotypes %>%
    replace(. == 0, -1) %>%
    replace(. == 2, 1) %>%
    replace(. == 3, 0) %>%
    as_tibble()

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
  "HDBSCAN Clusters", "Breeding Programs", "Era", "Market Class", "Phenotype",
  "Growth Habit", "Colour", "Texture"
)
categorization_colours <- list(
  colours_hdbscan_pic_legend, colours_bp, colours_era, colours_mc, colours_pheno,
  colour_set[c(1, 22, 4)], colour_set[c(1, 22, 4)], colour_set[c(1, 4, 22)]
)

lapply(seq_along(categorizations), function (i) {
  categorization <- categorizations[[i]]
  print(categorization_names[i])
  categories <- categorization %>% as.factor() %>% levels()

  rarefied <- mclapply(seq_along(categories), function (j) {
    category <- categories[j]
    individuals <- which(categorization == category)
    subset_sizes <- lseq(2, length(individuals), 10) %>%
      floor() %>% unique()
    # subset_sizes <- lseq(2, min(60, length(individuals)), 15) %>%
    #   floor() %>% unique()
    
    lapply(subset_sizes, function (subset_size) {
      print(str_c(category, ": ", subset_size))
      temp <- tibble(
        category = character(), subset_size = double(), pic = double()
      )
      for (i in 1:1000) {
        temp <- temp %>% add_row(
          category = categories[j], subset_size = subset_size,
          pic = PicCalc(wheat_formatted[, sample(individuals, subset_size)])
        )
      }
      temp
    }) %>% do.call(rbind, .)
  }, mc.cores = detectCores()) %>% do.call(rbind, .)

  plot <- rarefied %>% ggplot() +
    geom_point(
      aes(log10(subset_size), pic, colour = str_wrap(category, 10)), size = 0.5, alpha = 0.2,
      pch = 16
    ) +
    geom_smooth(
      aes(log10(jitter(subset_size)), pic, colour = str_wrap(category, 10)),
      se = FALSE, alpha = 0.5
    ) +
    scale_colour_manual(values = categorization_colours[[i]]) +
    labs(colour = "Category") +
    ggtitle(
      str_c(
        "Rarefaction Curves of Average PIC of\nSub-samples By ",
        categorization_names[i]
      )
    ) +
    xlab("Sample Size") +
    ylab("Average PIC") +
    guides(
      colour = guide_legend(nrow = (length(categories) / 3) %>% ceiling())) +
    theme(
      legend.position = "bottom", 
      legend.text = element_text(
        size = ifelse(length(categories) <= 4, 12,
          ifelse(length(categories) <= 8, 8, 4)
        )
      )
    )

  ggsave(
    str_c(
      "PIC_rarefaction_by_", str_replace(categorization_names[i], " ", "_"),
      ".png"
    ), plot = plot, path = PIC, width = 100, height = 120, units = "mm"
  )
})