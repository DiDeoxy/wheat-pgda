source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "arrange", "distinct", "mutate", "rowwise", "select")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2,
  "aes", "element_text", "ggplot", "geom_hline", "geom_line", "geom_point",
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
import::from(stringr, "str_c", "str_wrap")
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
  arrange(pos_mb)

gene_ranges <- read_csv(file.path(intermediate, "gene_ranges.csv"))

genes_nearby_markers <- lapply(1:nrow(gene_ranges), function (row) {
  markers <- snp_data[
    which(
      snp_data$chrom == gene_ranges[row, ]$chrom &
      snp_data$pos_mb >= gene_ranges[row, ]$window_start &
      snp_data$pos_mb <= gene_ranges[row, ]$window_end
    ),
  ] %>%
    add_row(pos_mb = gene_ranges[row, ]$window_start) %>%
    add_row(pos_mb = gene_ranges[row, ]$window_end) %>% 
    arrange(pos_mb)
  list(
    markers = markers,
    group = gene_ranges[row, ]$genes
  )
})
genes_nearby_markers[[1]]

legend_title <- "MJAF, Jost's D,\nand Gene Type"
lables <- c(
  "CHRS MJAF", "CHRW MJAF", "CSWS MJAF", "Jost's D", "Phenotype Genes",
  "Resistance Genes"
)
lapply(genes_nearby_markers, function(gene_data) {
  file_name <- str_c(
    str_c(
      gene_data$markers$chrom[2],
      round(min(gene_data$markers$pos_mb), 0),
      round(max(gene_data$markers$pos_mb), 0),
      gene_data$group,
      sep = "_"
    ),
    ".png"
  )
  print(file_name)
  png(
    file.path(zoomed_marker_plots, file_name),
    family = "Times New Roman",
    width = 320, height = 240,
    units = "mm", res = 192
  )
  print(gene_data$markers %>%
    ggplot(aes(pos_mb, values)) +
    ylim(0, 1) +
    xlim(range(gene_data$markers$pos_mb)) +
    geom_point(aes(colour = value_type, shape = value_type), size = 3) +
    geom_hline(yintercept = top_quartile) +
    geom_text_repel(
      aes(
        label = ifelse(
          value_type == "Phenotype Genes" | value_type == "Resistance Genes",
          id, ""
        )
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
    labs(
      x = "Position in Mb", y = "MJAF and Jost's D",
      title = str_c(gene_data$markers$chrom[2])
      # title = str_c(gene_data$chrom[1], gene_data$group[1], sep = "_")
    )
  )
  dev.off()
})