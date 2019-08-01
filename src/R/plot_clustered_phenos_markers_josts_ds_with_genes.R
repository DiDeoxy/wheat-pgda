# load file paths, colours, and functions
source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "arrange", "bind_cols", "mutate", "rowwise", "select")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_point", "geom_segment", "guide_legend",
  "scale_colour_manual", "scale_linetype", "theme", "xlim", "ylim"
)
import::from(ggrepel, "geom_label_repel")
import::from(magrittr, "%>%", "%<>%")
import::from(pgda, "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(plyr, "rbind.fill")
import::from(readr, "read_rds", "read_csv", "write_csv")
import::from(Rfast, "rowMaxs")
import::from(scrime, "rowTables")
import::from(stringr, "str_c")
import::from(tibble, "add_column", "add_row", "as_tibble", "tibble")
import::from(tidyr, "gather", "spread")

# load the data from the gds object
wheat_data <- snpgds_parse(phys_gds)

# recode the genotypes and count
genos <- replace(wheat_data$genotypes, wheat_data$geno == 3, NA)

# load other data
josts_d <- read_rds(josts_d)

# set the genotype coding
coding <- c(0, 2)

################################################################################
# get the cluster mjafs

cluster <- read_rds(hdbscan)$cluster
cluster_genos <- list(
  `CHRS MJAF` = genos[,
    which(wheat_data$sample$annot$pheno == "HRS" & cluster == 5)
  ],
  `CHRW MJAF` = genos[,
    which(wheat_data$sample$annot$pheno == "HRW" & cluster == 1)
  ],
  `CSWS MJAF` = genos[,
    which(wheat_data$sample$annot$pheno == "SWS" & cluster == 2)
  ]
)

# get the group MJAFs
mja <- rowMaxs(rowTables(genos, coding))
mjafs_by_pop <- lapply(cluster_genos, function (sub_genos) {
  genos_counts <- rowTables(sub_genos, coding)
  max_genos <- genos_counts[cbind(seq_along(mja), mja)]
  max_genos / rowSums(genos_counts)
}) %>% bind_cols() %>% round(4)

################################################################################
# make the stat data dataframe

# gen and phys tibbles of the data so far
snp_data <- tibble(
  chrom = wheat_data$snp$chrom,
  id = wheat_data$snp$id,
  pos_mb = wheat_data$snp$pos / 1e6,
  josts_d = josts_d,
)

snp_data %<>%
  bind_cols(mjafs_by_pop)

################################################################################
# classify josts d values

snp_data %<>%
  rowwise() %>%
  mutate(class =
    c("CHRSD Jost's D", "CHRWD Jost's D", "CSWSD Jost's D")[
      which.max(
        c(
          sum(
            abs(`CHRS MJAF` - `CHRW MJAF`), abs(`CHRS MJAF` - `CSWS MJAF`)
          ),
          sum(
            abs(`CHRW MJAF` - `CHRS MJAF`), abs(`CHRW MJAF` - `CSWS MJAF`)
          ),
          sum(
            abs(`CSWS MJAF` - `CHRS MJAF`), abs(`CSWS MJAF` - `CHRW MJAF`)
          )
        )
      )
    ]
  )

################################################################################
# create a tibble of the genes and centromeres

landmarks <- rbind.fill(
  load_genes(
    file.path(blast, "selected_pheno.csv")
  ) %>% mutate(
    type = "Gene", base = 0.5, pos_mb = pos / 1e6
  ) %>%
    select(-pos),
  load_genes(
    file.path(blast, "selected_resi.csv")
  ) %>% mutate(
    type = "Gene", base = 0.5, pos_mb = pos / 1e6
  ) %>%
    select(-pos),
  cbind(
    id = "Centromere", type = "Centromere", base = 0.5,
    read_csv(file.path(intermediate, "centromeres.csv"))
  )
)

################################################################################
# print out some stats

# overall freqs
summary(snp_data$josts_d)
sum(snp_data$josts_d > mean(snp_data$josts_d)) / length(snp_data$josts_d)

# genome freqs
for (genome in c("A", "B", "D")) {
  median(snp_data$josts_d[which(grepl(genome, snp_data$chrom))]) %>%
    round(4) %>% str_c(genome, " median = ", .) %>% print()
}

# chromsome group freqs
for (chr_group in 1:7) {
  median(snp_data$josts_d[which(grepl(chr_group, snp_data$chrom))]) %>%
    round(4) %>% str_c("Chr ", chr_group, " median = ", .) %>% print()
}

# group freq and median values
for (class in unique(snp_data$class)) {
  (sum(snp_data$class == class) / length(snp_data$class) * 100) %>%
    round(2) %>% str_c("Class ", class, " percent = ", .) %>% print()
  median(snp_data$josts_d[which(snp_data$class == class)]) %>%
    round(2) %>% str_c("Class ", class, " median = ", .) %>% print()
}

################################################################################
# reorganise the data

snp_data %<>%
  select(chrom, id, pos_mb, josts_d, class) %>%
  split(.$chrom)

################################################################################
# calc the lengths of the different genomes and homoeologous sets

max_lengths <- span_by_chrom(
  wheat_data$snp$chrom, wheat_data$snp$pos
) %>% max_lengths() / 1e6

################################################################################

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
legend_title_1 <- "Marker Types"
legend_title_2 <- "Landmarks"
plots <- lapply(snp_data, function(chrom_data) {
  chrom <- chrom_data$chrom[1]

  chrom_landmarks <- landmarks[which(landmarks$chrom == chrom), ]

  ggplot() +
    ylim(0, 1) +
    xlim(
      0,
      max_lengths[[
        ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
      ]]
    ) +
    geom_point(
      aes(chrom_data$pos_mb, chrom_data$josts_d, colour = chrom_data$class),
      size = 1
    ) +
    geom_segment(
      aes(
        chrom_landmarks$pos_mb, pmin(chrom_landmarks$base, 0),
        xend = chrom_landmarks$pos_mb, yend = 1, linetype = chrom_landmarks$type
      ),
      size = 1
    ) +
    geom_label_repel(
      aes(
        chrom_landmarks$pos_mb, chrom_landmarks$base, label = chrom_landmarks$id
      ),
      size = 1.5, seed = 101
    ) +
    scale_colour_manual(
      legend_title_1, values = colours_comps_genes
    ) +
    scale_linetype(
      legend_title_2
    )
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = str_c(
    "Normalized Jost's D Value"
  ),
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Marker-by-Marker Jost's D values Coloured by Distinct Group",
  legend = c(1, 1)
)

# plot the matrix
png(
  file.path("results", "clustered_phenos_markers_josts_ds_with_genes.png"),
  family = "Times New Roman", width = 500, height = 250, pointsize = 5,
  units = "mm", res = 192
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()
