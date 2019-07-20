source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "bind_cols", "mutate", "rowwise", "select")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_point", "geom_vline", "labs", 
  "scale_colour_manual",
  "theme", "unit", "xlim", "ylim"
)
import::from(magrittr, "%>%", "%<>%")
import::from(pgda, "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(plyr, "rbind.fill")
import::from(readr, "read_csv", "read_rds", "write_csv")
import::from(Rfast, "rowMaxs")
import::from(scrime, "rowTables")
import::from(stringr, "str_c")
import::from(tibble, "add_column", "as_tibble", "tibble")
import::from(tidyr, "gather", "spread")

# load the data from the gds object
phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

# recode the genotypes and count
genos <- replace(phys_data$genotypes, phys_data$geno == 3, NA)

# set the genotype coding
coding <- c(0, 2)

# load some data
josts_d <- read_rds(josts_d)
allele_counts <- rowTables(genos, coding)
mjafs <- rowMaxs(allele_counts, value = TRUE) / rowSums(allele_counts)
phys_to_gen_order <- match(gen_data$snp$id, phys_data$snp$id)
window_ranges <- read_csv(file.path(intermediate, "windows.csv"))

# gen and phys tibbles of the data so far
stat_data <- list(
  phys = tibble(
    chrom = phys_data$snp$chrom,
    id = phys_data$snp$id,
    pos_mb = phys_data$snp$pos / 1e6,
    josts_d = josts_d,
    mjafs = mjafs
  ),
  gen = tibble(
    chrom = gen_data$snp$chrom,
    id = gen_data$snp$id,
    pos_cm = gen_data$snp$pos / 100,
    josts_d = josts_d[phys_to_gen_order],
    mjafs = mjafs[phys_to_gen_order]
  )
)

################################################################################
# caclualte each markers order difference between gen and phys maps and assign
# it to an interval

order_diffs <- by(
  bind_cols(phys_id = stat_data$phys$id, gen_id = stat_data$gen$id),
  stat_data$phys$chrom,
  function (chrom_data) {
    # find the relative order of gen markers to phys markers
    gen_to_phys_order <- match(chrom_data$phys_id, chrom_data$gen_id)

    # find the magnitude of the order differences and turn into intervals
    (gen_to_phys_order - 1:length(gen_to_phys_order)) %>% abs()
  }
) %>% unlist()

order_diff_quantiles <- c("0%" = -1,
  quantile(
    order_diffs, c(1 / 2, 3 / 4, 7 / 8, 15 / 16, 31 / 32, 63 / 64, 127 / 128)
  ),
  "100%" = max(order_diffs)
)

intervals <- lapply(
  seq_along(order_diff_quantiles),
  function (i) {
    if (i < length(order_diff_quantiles)) {
      str_c(order_diff_quantiles[i] + 1, "-", order_diff_quantiles[i + 1])
    }
  }
) %>% unlist()

order_diff_intervals <- cut(order_diffs, order_diff_quantiles, intervals)

stat_data$phys %<>%
  add_column(type = order_diff_intervals)
stat_data$gen %<>%
  add_column(type = order_diff_intervals[phys_to_gen_order])



################################################################################
# get the cluster mjafs

cluster <- read_rds(hdbscan)$cluster
cluster_genos <- list(
  `CHRS MJAF` = genos[,
    which(phys_data$sample$annot$pheno == "HRS" & cluster == 5)
  ],
  `CHRW MJAF` = genos[,
    which(phys_data$sample$annot$pheno == "HRW" & cluster == 1)
  ],
  `CSWS MJAF` = genos[,
    which(phys_data$sample$annot$pheno == "SWS" & cluster == 2)
  ]
)

# get the group MJAFs
mja <- rowMaxs(rowTables(genos, coding))
mjafs_by_pop <- lapply(cluster_genos, function (sub_genos) {
  genos_counts <- rowTables(sub_genos, coding)
  max_genos <- genos_counts[cbind(seq_along(mja), mja)]
  max_genos / rowSums(genos_counts)
}) %>% bind_cols()

stat_data$phys %<>%
  bind_cols(mjafs_by_pop %>% round(4))

stat_data$gen %<>%
  bind_cols(mjafs_by_pop[phys_to_gen_order, ] %>% round(4))

################################################################################
# add in genes

all_genes <- rbind.fill(
  load_genes(
    file.path(blast, "selected_pheno.csv")
  ) %>% mutate(type = "Phenotype Gene", base = 1, pos_mb = pos / 1e6) %>%
    select(-pos),
  load_genes(
    file.path(blast, "selected_resi.csv")
  ) %>% mutate(type = "Resistance Gene", base = 1, pos_mb = pos / 1e6) %>%
    select(-pos)
)

stat_data$phys %<>%
  rbind.fill(all_genes)

stat_data$gen %<>%
  rbind.fill(all_genes)

################################################################################
# classify josts d values

stat_data$phys %<>%
    rowwise() %>%
    mutate(class =
      c("CHRSD Jost's D", "CHRWD Jost's D", "CSWSD Jost's D", "None")[
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
            ),
            0
          )
        )
      ]
    )
# by(
#   stat_data$phys, stat_data$phys$chrom, function (chrom_data) {
#     write_csv(
#       chrom_data %>% select(
#         id, chrom, pos_mb, josts_d, `CHRS MJAF`, `CHRW MJAF`, `CSWS MJAF`, class
#       ),
#       file.path(josts_d_by_chrom,
#         str_c(
#           chrom_data$chrom[1],
#           "_all_markers_josts_d_and_mjaf_freqs_by_pop_with_genes.csv"
#         )
#       )
#     )
#   }
# )
stat_data$phys %<>% spread(class, josts_d)

stat_data$gen %<>%
    rowwise() %>%
    mutate(class =
      c("CHRSD Jost's D", "CHRWD Jost's D", "CSWSD Jost's D", "None")[
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
            ),
            0
          )
        )
      ]
    ) %>%
    spread(class, josts_d)

################################################################################
# calc each groups mjaf

# add the centromere positions
stat_data$phys %<>%
  rbind.fill(
    cbind(
      type = "Centromere",
      read_csv(
        file.path(intermediate, "centromeres.csv")
      )[, c(1, 2)]
    )
  )

stat_data$gen %<>%
  rbind.fill(
    cbind(
      type = "Centromere",
      read_csv(
        file.path(intermediate, "centromeres.csv")
      )[, c(1, 3)]
    )
  )

################################################################################
# reorganize the data

stat_data$phys %<>%
  mutate(type = as.factor(type)) %>%
  gather(stat, value, -c(chrom, id, pos_mb, cent_pos_mb, mjafs, type, base)) %>%
  split(.$chrom)

stat_data$gen %<>%
  mutate(type = as.factor(type)) %>%
  gather(stat, value, -c(chrom, id, pos_mb, cent_pos_cm, mjafs, type, base)) %>%
  split(.$chrom)

head(stat_data$phys[[1]])
unique(stat_data$phys[[1]]$type)
unique(stat_data$phys[[1]]$stat)

stat_data$phys[[1]][which(stat_data$phys[[1]]$type == "Resistance Gene"), ]

################################################################################
# plot the data

lapply(stat_data$phys[1], function (chrom_data) {
  chrom <- chrom_data$chrom[1]
  chrom_window_ranges <- window_ranges[which(windows_data$chrom == chrom), ]

  mjafs <- chrom_data %>% ggplot() +
    geom_point(
      aes(
        pos_mb, mjafs, colour = type
      ),
      size = 2
    ) +
    geom_vline(
      aes(xintercept = cent_pos_mb, colour = type), size = 2
    ) +
    scale_colour_manual(
      name = "Order Difference\n& Centromere",
      values = colours_intervals,
      drop = FALSE
    )

  cluster_mjafs_josts_d_1 <- chrom_data %>% ggplot(aes(pos_mb, value)) +
    ylim(0, 1) +
    geom_point(
      aes(
        colour = type, shape = type,
        size = type, fill = type
      )
    ) +
    geom_text_repel(
      aes(
        label = ifelse(
          type == "Phenotype Genes" | type == "Resistance Genes",
          id, NA
        )
      ),
      size = 5
    ) +
    scale_fill_manual(
      legend_title, #labels = lables,
      values = colour_set[c(1, 15, 2, 17, 4, 19, 15, 19)],
      na.translate = FALSE
    ) +
    scale_shape_manual(
      legend_title, #labels = lables,
      values = c(16, 22, 16, 22, 16, 22, 25, 25),
      na.translate = FALSE
    ) +
    scale_size_manual(
      legend_title, #labels = lables,
      values = c(2, 3, 2, 3, 2, 3, 4, 4),
      na.translate = FALSE
    ) +
    scale_colour_manual(
      legend_title, #labels = lables,
      values = colour_set[c(1, 22, 2, 22, 4, 22, 15, 19)],
      na.translate = FALSE
    ) +
    scale_x_continuous(breaks = seq(0, max(chrom_data$pos_mb), by = 10)) +
    labs(
      x = "Position in Mb", y = "MJAF and Jost's D",
      title = str_c(chrom_data$chrom[1])
    )
  cluster_mjafs_josts_d_2 <- NA
  if (nrow(chrom_window_ranges)) {
    for (row in 1:nrow(chrom_window_ranges)) {
      if (any(is.na(cluster_mjafs_josts_d_2))) {
        cluster_mjafs_josts_d_2 <- cluster_mjafs_josts_d +
          geom_rect(
            data = chrom_window_ranges[row, ], inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = colour_set[20], alpha = 0.3
          )
      } else {
        cluster_mjafs_josts_d_2 <- cluster_mjafs_josts_d_2 +
          geom_rect(
            data = chrom_window_ranges[row, ], inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = colour_set[20], alpha = 0.3
          )
      }
      # zoomed plots
      file_name <- str_c(
        str_c(
          chrom,
          round(chrom_window_ranges[row, ]$start, 0),
          round(chrom_window_ranges[row, ]$end, 0),
          chrom_window_ranges[row, ]$genes,
          sep = "_"
        ),
        ".png"
      )
      print(file_name)
      png(
        file.path(zoomed_marker_plots, file_name),
        family = "Times New Roman",
        width = 300, height = 150,
        units = "mm", res = 192
      )
      print(
        cluster_mjafs_josts_d_1 +
          xlim(chrom_window_ranges[row, ]$start, chrom_window_ranges[row, ]$end)
      )
      dev.off()
    }
  }

  chrom_data %>% ggplot() +


  if (! any(is.na(cluster_mjafs_josts_d_2))) {
    return(
      list(
        mjafs, cluster_mjafs_josts_d_2
      )
    )
  } else {
    return(
      list(
        mjafs, cluster_mjafs_josts_d_1
      )
    )
  }
})

png(
  file.path(all_data_chroms, str_c(chrom, ".png")),
  family = "Times New Roman",
  width = 1800, height = 300,
  units = "mm", res = 192
)
dev.off()