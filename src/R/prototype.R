source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "arrange", "distinct", "mutate", "rowwise", "select")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2,
  "aes", "element_text", "ggplot", "geom_hline", "geom_line", "geom_rect",
  "geom_point",
  "geom_smooth",
  "guide_legend", "labs", "scale_colour_manual", "scale_shape_manual", 
  "scale_x_continuous", "theme", "xlab", "ylab", "xlim", "ylim"
)
import::from(ggrepel, "geom_text_repel")
import::from(magrittr, "%>%")
import::from(
  pgda, "calc_eh", "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom"
)
import::from(plyr, "rbind.fill")
import::from(readr, "read_csv", "read_rds", "write_csv")
import::from(Rfast, "rowMaxs")
import::from(scrime, "rowTables")
import::from(stringr, "str_c", "str_replace", "str_wrap")
import::from(tibble, "add_column", "add_row", "as_tibble", "tibble")
import::from(tidyr, "gather", "spread")

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
  mutate(pos_mb = pos / 1e6) %>%
  arrange(chrom, pos_mb) %>%
  select(-pos) %>%
  add_column(josts_d := read_rds(josts_d), class = "josts_d") %>%
  cbind(mjafs_by_pop)

top_quartile <- snp_data$josts_d %>% quantile(0.75, na.rm = T)

# load the gene positions
all_genes <- rbind.fill(
  load_genes(
    file.path(blast, "selected_pheno.csv"), base = 1
  ) %>% mutate(gene_type = "Phenotype Genes", pos_mb = pos / 1e6) %>%
    select(-pos),
  load_genes(
    file.path(blast, "selected_resi.csv"), base = 1
  ) %>% mutate(gene_type = "Resistance Genes", pos_mb = pos / 1e6) %>%
    select(-pos)
) %>% arrange(chrom, pos_mb)

snp_data <- snp_data %>%
  rbind.fill(all_genes) %>%
  mutate(
    type = pmin(gene_type, class, na.rm = TRUE),
    jdb = pmin(josts_d, base, na.rm = TRUE)
  ) %>%
  select(-c(gene_type, class, base, josts_d)) %>%
  spread(type, jdb) %>%
  gather(value_type, values, -c(id, chrom, pos_mb)) %>%
  arrange(chrom, pos_mb)

chroms_data <- split(snp_data, snp_data$chrom)
windows_data <- read_csv(file.path(intermediate, "windows.csv"))

legend_title <- "MJAF, Jost's D,\nand Gene Type"
lables <- c(
  "CHRS MJAF", "CHRW MJAF", "CSWS MJAF", "Jost's D", "Phenotype Genes",
  "Resistance Genes"
)
lapply(chroms_data[10], function(chrom_data) {
  chrom <- chrom_data$chrom[1]
  print(chrom)
  chrom_windows <- windows_data[which(windows_data$chrom == chrom), ]

  # full chrom plot
  p1 <- chrom_data %>%
    ggplot(aes(pos_mb, values)) +
    ylim(0, 1) +
    geom_point(
      aes(colour = value_type, shape = value_type),
      size = 3
    ) +
    geom_hline(yintercept = top_quartile) +
    geom_text_repel(
      aes(
        label = ifelse(
          value_type == "Phenotype Genes" | value_type == "Resistance Genes",
          id, NA
        ),
        label_size = 1, xlim = range(chrom_data$pos_mb)
      )
    ) +
    scale_colour_manual(
      legend_title, labels = lables,
      values = colour_set[c(1, 2, 4, 22, 15, 19)],
      na.translate = FALSE
    ) +
    scale_shape_manual(
      legend_title, labels = lables,
      values = c(16, 16, 16, 15, 25, 25),
      na.translate = FALSE
    ) +
    scale_x_continuous(breaks = seq(0, max(chrom_data$pos_mb), by = 10)) +
    labs(
      x = "Position in Mb", y = "MJAF and Jost's D",
      title = str_c(chrom_data$chrom[1])
      # title = str_c(chrom_data$chrom[1], chrom_data$group[1], sep = "_")
    )
  p2 <- NA
  if (nrow(chrom_windows)) {
    for (row in 1:nrow(chrom_windows)) {
      if (any(is.na(p2))) {
        p2 <- p1 +
        geom_rect(
          data = chrom_windows[row, ], inherit.aes = FALSE,
          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
          fill = colour_set[20], alpha = 0.3
        )
      } else {
        p2 <- p2 +
        geom_rect(
          data = chrom_windows[row, ], inherit.aes = FALSE,
          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
          fill = colour_set[20], alpha = 0.3
        )
      }

      # zoomed plots
      file_name <- str_c(
        str_c(
          chrom,
          round(chrom_windows[row, ]$start, 0),
          round(chrom_windows[row, ]$end, 0),
          chrom_windows[row, ]$genes,
          sep = "_"
        ),
        ".png"
      )
      print(file_name)
      png(
        file.path(zoomed_marker_plots, file_name),
        family = "Times New Roman",
        width = 600, height = 200,
        units = "mm", res = 192
      )
      print(
        p1 +
          xlim(chrom_windows[row, ]$start, chrom_windows[row, ]$end)
      )
      dev.off()
    }
  }
  png(
    file.path(chrom_mjaf_josts_d, str_c(chrom_data$chrom[1], ".png")),
    family = "Times New Roman",
    width = 1800, height = 300,
    units = "mm", res = 192
  )
  print(p2)
  dev.off()
})