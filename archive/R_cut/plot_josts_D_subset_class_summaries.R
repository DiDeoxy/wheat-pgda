source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "arrange", "mutate", "rowwise", "select")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2,
  "aes", "element_text", "ggplot", "geom_line", "geom_point",
  "guide_legend", "labs", "scale_colour_manual", "scale_x_continuous", "theme",
  "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(pgda, "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(plyr, "rbind.fill")
import::from(readr, "read_rds", "write_csv")
import::from(Rfast, "rowMaxs")
import::from(scrime, "rowTables")
import::from(stringr, "str_c", "str_wrap")
import::from(tibble, "add_column", "add_row", "as_tibble", "tibble")
import::from(tidyr, "gather")

# load the data from the gds object
wheat_data <- snpgds_parse(phys_gds)

# recode the genotypes
wheat_data$genotypes <- replace(wheat_data$genotypes, wheat_data$geno == 3, NA)

# get the cluster groups
cluster <- read_rds(hdbscan)$cluster
genos <- list(
  chrs = wheat_data$genotype[,
    which(wheat_data$sample$annot$pheno == "HRS" & cluster == 5)
  ],
  chrw = wheat_data$genotype[,
    which(wheat_data$sample$annot$pheno == "HRW" & cluster == 1)
  ],
  csws = wheat_data$genotype[,
    which(wheat_data$sample$annot$pheno == "SWS" & cluster == 2)
  ]
)

# calc each markers mjaf by each cluster group
coding <- c(0, 2)
mja <- rowMaxs(rowTables(wheat_data$genotypes, coding))
mjafs_by_pop <- lapply(genos, function (geno) {
  geno_counts <- rowTables(geno, coding)
  max_genos <- geno_counts[cbind(seq_along(mja), mja)]
  max_genos / rowSums(geno_counts)
}) %>% do.call(cbind, .)

# add the genes positons to the regions table
snp_data <- wheat_data$snp %>%
  add_column(josts_d := read_rds(josts_d)) %>%
  cbind(mjafs_by_pop) %>% as_tibble()
top_quartile <- snp_data$josts_d %>% quantile(0.75, na.rm = T)
snp_data <- snp_data %>%
  rowwise() %>%
  mutate(class =
    ifelse(
      josts_d > top_quartile,
      c("CHRSD", "CHRWD", "CSWSD")[
        which.max(
          c(
            sum(abs(chrs - chrw), abs(chrs - csws)),
            sum(abs(chrw - chrs), abs(chrw - csws)),
            sum(abs(csws - chrs), abs(csws - chrw))
          )
        )
      ],
      "None"
    )
  )

# calc and print out freqs
all_D <- summary(snp_data$josts_d) %>% as.vector()
summaries <- tibble(
  "Subset" = "All",
  "Min." = all_D[1], "1st Qu." = all_D[2], "Median" = all_D[3],
  "Mean" = all_D[4], "3rd Qu." = all_D[5], "Max." = all_D[6],
  "Fraction of Markers" = 1
)
for (genome in c("A", "B", "D", 1:7)) {
  markers_on_genome <- which(grepl(genome, snp_data$chrom))
  genome_D <- summary(snp_data$josts_d[markers_on_genome])

  summaries <- summaries %>%
    add_row(
      Subset = genome,
      "Min." = genome_D[1], "1st Qu." = genome_D[2], "Median" = genome_D[3],
      "Mean" = genome_D[4], "3rd Qu." = genome_D[5], "Max." = genome_D[6],
      "Fraction of Markers" = length(markers_on_genome) / nrow(snp_data)
    )
}

# group freq and median values
classes <- c("CHRSD", "CHRWD", "CSWSD", "None")

for (class in classes) {
  markers_in_class <- which(snp_data$class == class)
  class_D <- summary(snp_data$josts_d[markers_in_class])

  summaries <- summaries %>%
    add_row(
      Subset = class,
      "Min." = class_D[1], "1st Qu." = class_D[2], "Median" = class_D[3],
      "Mean" = class_D[4], "3rd Qu." = class_D[5], "Max." = class_D[6],
      "Fraction of Markers" = length(markers_in_class) / nrow(snp_data)
    )
}
summaries[, 2:8] <- summaries[, 2:8] %>% round(3)

# print out summaries
write_csv(summaries, file.path("results", "josts_D_subset_class_summaries.csv"))

summaries_long <- summaries %>% gather("Statistic", "Values", -Subset)
summaries_long$Statistic <- factor(
  summaries_long$Statistic, levels = unique(summaries_long$Statistic)
)
summaries_long$Subset <- factor(
  summaries_long$Subset, levels = rev(unique(summaries_long$Subset))
)

# graph summaries
png(
  file.path("results", "josts_D_subset_class_summaries.png"),
  family = "Times New Roman", width = 208, height = 104, units = "mm", res = 300
)
ggplot(summaries_long) +
  geom_point(aes(Values, Subset, colour = Statistic), size = 1.5) +
  scale_colour_manual(
    values = colour_set[c(1, 5, 3, 2, 4, 6, 22)],
    labels = str_wrap(levels(summaries_long$Statistic), 10)
  ) +
  scale_x_continuous(breaks = 0:10 / 10)
dev.off()