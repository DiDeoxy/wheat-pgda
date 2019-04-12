source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "bind_rows")
import::from(emdbook, "lseq")
import::from(GeneticSubsetter, "PicCalc")
import::from(
  ggplot2, "aes", "geom_point", "geom_smooth", "ggplot", "ggtitle", "guide_legend",
  "guides", "labs", "scale_colour_manual", "theme", "xlab", "ylab"
)
import::from(magrittr, "%>%")
import::from(pgda, "snpgds_parse")
import::from(readr, "read_rds")
import::from(stringr, "str_c")
import::from(tibble, "add_row", "tibble")

wheat_data <- snpgds_parse(phys_gds)

wheat_formatted <- wheat_data$genotypes %>%
    replace(. == 0, -1) %>%
    replace(. == 2, 1) %>%
    replace(. == 3, 0)

clusters <- factor(read_rds(hdbscan)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1 (HRW)", "Cluster 2 (SWS)", "Cluster 3 (CWES)",
  "Cluster 4 (CPSR/W)", "Cluster 5 (HRS)"
)
levels(clusters)

rarefied <- tibble(cluster = character(), subset_size = double(), pic = double())
for (cluster in levels(clusters)[2:length(clusters)]) {
  print(cluster)
  samples <- which(clusters == cluster)
  subset_sizes <- lseq(2, length(samples) - 1, 20) %>% floor() %>% unique()
  for (subset_size in subset_sizes) {
    print(subset_size)
    for (i in 1:100) {
      subset <- sample(samples, subset_size)
      rarefied <- rarefied %>% 
        add_row(
          cluster = cluster, subset_size = subset_size, 
          pic = PicCalc(wheat_formatted[, subset])
        )
    }
  }
}

png(
  file.path("results", "PIC_rarefaction_by_cluster.png"),
  family = "Times New Roman", width = 100, height = 120, pointsize = 5,
  units = "mm", res = 300)
rarefied %>% ggplot() +
  geom_point(aes(log10(subset_size), pic, colour = cluster), size = 0.5, alpha = 0.4, pch = 16) +
  geom_smooth(aes(log10(jitter(subset_size)), pic, colour = cluster), se = FALSE) +
  scale_colour_manual(values = colours_hdbscan_legend[2:length(colours_hdbscan_legend)]) +
  labs(colour = "Cluster") +
  ggtitle("Rarefaction Curves of Average PIC of\nSub-samples By Cluster") +
  xlab("Log 10 Sample Size") +
  ylab("Average PIC") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))
dev.off()

groups <- str_c(wheat_data$sample$annot$era, wheat_data$sample$annot$bp)
groups <- str_c(wheat_data$sample$annot$pheno)
rarefied <- tibble(group = character(), subset_size = double(), pic = double())
for (group in unique(groups)) {
  print(group)
  samples <- which(groups == group)
  if (length(samples) >= 5) {
    subset_sizes <- lseq(2, length(samples) - 1, 20) %>% floor() %>% unique()
    for (subset_size in subset_sizes) {
      print(subset_size)
      for (i in 1:100) {
        subset <- sample(samples, subset_size)
        rarefied <- rarefied %>% 
          add_row(
            group = group, subset_size = subset_size, 
            pic = PicCalc(wheat_formatted[, subset])
          )
      }
    }
  }
}

rarefied %>% ggplot() +
  geom_point(aes(log10(subset_size), pic, colour = group), size = 0.5, alpha = 0.4, pch = 16) +
  geom_smooth(aes(log10(jitter(subset_size)), pic, colour = group), se = FALSE, n = 100, span = 1) +
  # scale_colour_manual(values = colours_hdbscan_legend[2:length(colours_hdbscan_legend)]) +
  labs(colour = "group") +
  ggtitle("Rarefaction Curves of Average PIC of\nSub-samples By group") +
  xlab("Log 10 Sample Size") +
  ylab("Average PIC") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 5))