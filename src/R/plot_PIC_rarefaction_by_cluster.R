source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "bind_rows")
import::from(emdbook, "lseq")
import::from(GeneticSubsetter, "PicCalc")
import::from(
  ggplot2, "aes", "annotation_logticks", "element_blank", "element_text",
  "geom_point", "geom_line", "ggplot", "ggsave", "ggtitle","guide_legend",
  "guides", "labs", "scale_colour_manual", "scale_x_continuous", "theme",
  "xlab", "ylab"
)
import::from(magrittr, "%>%")
import::from(parallel, "detectCores", "mclapply")
import::from(pgda, "snpgds_parse")
import::from(readr, "read_rds")
import::from(stringr, "str_c", "str_replace", "str_wrap")
import::from(tibble, "as_tibble", "tibble")

wheat_data <- snpgds_parse(phys_gds)

wheat_formatted <- wheat_data$genotypes %>%
    replace(. == 0, -1) %>%
    replace(. == 2, 1) %>%
    replace(. == 3, 0) %>%
    as_tibble()

clusters <- factor(read_rds(hdbscan)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1 (HRW)", "Cluster 2 (SWS)", "Cluster 3 (CWES)",
  "Cluster 4 (CPSR/W)", "Cluster 5 (CXRS)"
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
  colours_hdbscan_pic_legend, colours_bp_pic, colours_era, colours_mc, colours_pheno_pic,
  colour_set[c(1, 22, 4)], colour_set[c(1, 22, 4)], colour_set[c(1, 4, 22)]
)

lapply(4:8, function (i) {
  categorization <- categorizations[[i]]
  print(categorization_names[i])
  categories <- categorization %>% as.factor() %>% levels()

  rarefied <- mclapply(seq_along(categories), function (j) {
    category <- categories[j]
    individuals <- which(categorization == category)
    subset_sizes <- lseq(2, max(2, length(individuals)), 20) %>%
      floor() %>% unique()
    if (length(subset_sizes) > 1) {
      lapply(subset_sizes, function (subset_size) {
        print(str_c(category, ": ", subset_size))
        tibble(
          category = categories[j], subset_size = subset_size,
          pic = mean(
            sapply(1:1000, function (k) {
              PicCalc(wheat_formatted[, sample(individuals, subset_size)])
            })
          )
        )
      }) %>% do.call(rbind, .)
    }
  }, mc.cores = detectCores()) %>% do.call(rbind, .)

  plot <- rarefied %>% ggplot() +
    geom_point(
      aes(log10(subset_size), pic, colour = str_wrap(category, 10)), pch = 16
    ) +
    geom_line(
      aes(subset_size %>% log10(), pic, colour = str_wrap(category, 10))
    ) +
    scale_x_continuous(
      labels = scales::math_format(10^.x)
    ) +
    scale_colour_manual(values = categorization_colours[[i]]) +
    ggtitle(
      str_c(
        "Rarefaction Curves of Average PIC of\nSub-samples By ",
        categorization_names[i]
      )
    ) +
    xlab("Sample Size") +
    ylab("Average PIC") +
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
    annotation_logticks()

  ggsave(
    str_c(
      "PIC_rarefaction_by_", str_replace(categorization_names[i], " ", "_"),
      ".png"
    ), plot = plot, path = PIC, width = 100, height = 120, units = "mm"
  )
})