source("wheat-pgda/src/R/file_paths.R")
source("wheat-pgda/src/R/colour_sets.R")
import::from(ape, "read.gff")
import::from(
  dplyr, "arrange", "bind_cols", "filter", "left_join", "mutate", "mutate_at",
  "rename", "rowwise", "select", "ungroup"
)
import::from(GGally, "ggmatrix")
library(ggplot2)
import::from(ggpubr, "as_ggplot", "get_legend")
import::from(ggrepel, "geom_label_repel")
import::from(gridExtra, "grid.arrange")
import::from(magrittr, "%>%", "%<>%")
import::from(pgda, "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(plyr, "rbind.fill")
import::from(
  readr, "col_character", "col_double", "col_factor", "col_integer", "read_csv", 
  "read_rds", "type_convert", "write_csv"
)
import::from(scrime, "rowTables")
import::from(stringr, "str_c")
import::from(tibble, "add_column", "as_tibble", "enframe", "tibble")
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
mjafs <- do.call(pmax, (allele_counts / rowSums(allele_counts)) %>% as.data.frame())
eh <- (mjafs * (1 - mjafs)) * 2
gen_to_phys_order <- match(phys_data$snp$id, gen_data$snp$id)
chroms <- as.vector(
  t(outer(as.character(1:7), c("A", "B", "D"), paste, sep = ""))
)

################################################################################
# caclualte each markers order difference between gen and phys maps 

gen_to_phys_order <- match(phys_data$snp$id, gen_data$snp$id)
order_diffs <- (gen_to_phys_order - 1:length(gen_to_phys_order)) %>% abs()

png(file.path("results", "test.png"))
order_diffs %>% enframe() %>% ggplot(aes(value)) +
  # stat_ecdf() +
  geom_density() +
  scale_x_continuous(
    trans = "pseudo_log",
    breaks = c(0, 1, 2, 3, 4, 5, 10, 25, 50, 100, 250, 500, 1000)
  ) +
  labs(y = "Cumulative Density", x = "Marker Order Difference")
dev.off()

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

################################################################################
# get the cluster mjafs

cluster <- read_rds(hdbscan)$cluster
cluster_genos <- list(
  CHRS = genos[,
    which(phys_data$sample$annot$mtg == "HRS" & cluster == 5)
  ],
  CHRW = genos[,
    which(phys_data$sample$annot$mtg == "HRW" & cluster == 1)
  ],
  CSWS = genos[,
    which(phys_data$sample$annot$mtg == "SWS" & cluster == 2)
  ]
)

# get the group MJAFs
mja <- apply(rowTables(genos, coding), 1, which.max)
cluster_mjafs <- lapply(cluster_genos, function (sub_genos) {
  genos_counts <- rowTables(sub_genos, coding)
  max_genos <- genos_counts[cbind(seq_along(mja), mja)]
  max_genos / rowSums(genos_counts)
}) %>% bind_cols()

################################################################################
# gen and phys tibbles of the data so far
snp_data <- tibble(
  chrom = phys_data$snp$chrom,
  id = phys_data$snp$id,
  pos_mb = phys_data$snp$pos / 1e6,
  pos_cm = gen_data$snp$pos[gen_to_phys_order] / 100,
  eh = eh,
  odi = factor(order_diff_intervals),
  josts_d = josts_d
) %>%
  bind_cols(cluster_mjafs) %>%
  rowwise() %>%
  mutate(
    josts_d_class = c("CHRSD", "CHRWD", "CSWSD")[
      which.max(
        c(
          sum(abs(CHRS - CHRW), abs(CHRS - CSWS)),
          sum(abs(CHRW - CHRS), abs(CHRW - CSWS)),
          sum(abs(CSWS - CHRS), abs(CSWS - CHRW))
        )
      )
    ]
  ) %>%
  ungroup() %>%
  mutate(josts_d_class = factor(josts_d_class)) %>%
  split(.$chrom)

# # mjafs and ehs by cluster
# mjafs_by_cluster <- cluster_mjafs %>%
#   bind_cols(chrom = phys_data$snp$chrom, pos_mb = phys_data$snp$pos / 1e6) %>%
#   gather(cluster, mjaf, CHRS, CHRW, CSWS) %>%
#   as_tibble() %>%
#   split(.$chrom)

# ehs_by_cluster <- ((cluster_mjafs * (1 - cluster_mjafs)) * 2) %>%
#   bind_cols(chrom = phys_data$snp$chrom, pos_mb = phys_data$snp$pos / 1e6) %>%
#   gather(cluster, eh, CHRS, CHRW, CSWS) %>%
#   as_tibble() %>%
#   split(.$chrom)

################################################################################
# create a tibble of the genes and centromeres

landmarks <- rbind.fill(
  load_genes(
    file.path(blast, "selected_pheno.csv")
  ) %>% mutate(
    type = "Gene", base = 0.75, pos_mb = pos / 1e6
  ) %>%
    dplyr::select(-pos),
  load_genes(
    file.path(blast, "selected_resi.csv")
  ) %>% mutate(
    type = "Gene", base = 0.75, pos_mb = pos / 1e6
  ) %>%
    dplyr::select(-pos),
  cbind(
    id = "Centromere", type = "Centromere", base = 0.75,
    read_csv(file.path(intermediate, "centromeres.csv"))
  )
) %>% as_tibble()

################################################################################
# or use go genes as landmarks

# landmarks <- read_csv(
#   file.path("/workspace", "data", "intermediate", "auxin_genes.csv")
# )

# landmarks <- read_csv(
#   file.path("/workspace", "data", "intermediate", "wounding_genes.csv")
# )

################################################################################
# calc the lengths of the different genomes and homoeologous sets

max_phys_lengths <- span_by_chrom(
  phys_data$snp$chrom, phys_data$snp$pos
) %>% max_lengths() / 1e6
max_gen_lengths <- span_by_chrom(
  gen_data$snp$chrom, gen_data$snp$pos
) %>% max_lengths() / 100

################################################################################
# calc the distirbution of markers in mono haplos and those not

mono_haplos <- cbind(
  id = "mono_haplo", type = "mono_haplo", base = 0.5,
  read_csv(file.path(intermediate, "haplo_windows.csv"))
) %>%
  split(.$chrom)

plots <- lapply(names(mono_haplos), function (chrom) {
  chrom_data <- snp_data[[chrom]] %>%
    rowwise() %>%
    mutate(
      type = c("All Other", "Long Haplotype")[
        which(
          c(
            pos_mb < mono_haplos[[chrom]]$start |
            pos_mb > mono_haplos[[chrom]]$end,
            pos_mb >= mono_haplos[[chrom]]$start &
            pos_mb <= mono_haplos[[chrom]]$end
          )
        )
      ]
    )

  chrom_data %>% ggplot(aes(type, eh)) +
    ylim(0, 0.5) +
    geom_violin(aes(fill = type)) +
    geom_boxplot(width = 0.2) +
    labs(y = "Diversity", title = chrom) +
    theme(legend.position = "none")
})

png(
  file.path("results", "mono_haplo_vs_distal_eh_dists.png"),
  family = "Times New Roman", width = 300, height = 600, pointsize = 5,
  units = "mm", res = 192
)
grid.arrange(
  grobs = plots, nrow = 3, ncol = 2
)
dev.off()

################################################################################
# plot the data
line_size <- 1
point_size <- 2
text_size <- 10
seed <- 101

remove_x_axis <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

remove_y_axis <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

plots <- sapply(chroms, function (chrom) {
  print(chrom)

  chrom_landmarks <- landmarks[which(landmarks$chrom == chrom), ]

  max_pos_cm <- max_gen_lengths[[
    ifelse(grepl("1", chrom), "one",
      ifelse(grepl("2", chrom), "two",
        ifelse(grepl("3", chrom), "three",
          ifelse(grepl("4", chrom), "four",
            ifelse(grepl("5", chrom), "five",
              ifelse(grepl("6", chrom), "six", "seven")
            )
          )
        )
      )
    )
  ]]

  max_pos_mb <- max_phys_lengths[[
    ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
  ]]

  # list for point_plots
  point_plots <- list()
  dist_plots <- list()

  plot_landmarks <- function(max_y) {
    list(
      geom_vline(
        aes(
          xintercept = chrom_landmarks$pos_mb,
          linetype = chrom_landmarks$type
        ), size = line_size
      ),
      geom_label_repel(
        aes(
          chrom_landmarks$pos_mb, chrom_landmarks$base * max_y,
          label = chrom_landmarks$id
        ), size = text_size, seed = seed, inherit.aes = FALSE
      ),
      scale_linetype("Landmarks")
    )
  }

  ## 1
  dist_plots$mjafs <- ggplot(mjafs_by_cluster[[chrom]], aes(cluster, mjaf)) +
    geom_violin(aes(fill = cluster)) +
    geom_boxplot(width = 0.2) +
    expand_limits(y = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_fill_manual(
      "Cluster MJAF", values = colour_set[c(1, 2, 4)],
      na.translate = FALSE, drop = FALSE
    ) +
    remove_x_axis +
    labs(y = "Major Allele Frequency", title = chrom) +
    guides(fill = FALSE)

  point_plots$mjafs <- ggplot() +
    scale_x_continuous(
      breaks = seq(0, max_pos_mb, by = 10), expand = c(0.01, 0.01)
    ) +
    expand_limits(x = c(0, max_pos_mb), y = c(0, 1)) +
    geom_point(
      aes(
        mjafs_by_cluster[[chrom]]$pos_mb, mjafs_by_cluster[[chrom]]$mjaf,
        colour = mjafs_by_cluster[[chrom]]$cluster
      ), size = point_size
    ) +
    scale_colour_manual(
      "Cluster MJAF",
      values = colour_set[c(1, 2, 4)], na.translate = FALSE, drop = FALSE
    ) +
    plot_landmarks(1) +
    labs(y = "Major Allele Frequency", title = chrom) +
    remove_x_axis +
    guides(linetype = FALSE)

  ## 2
  dist_plots$cluster_eh <- ggplot(ehs_by_cluster[[chrom]], aes(cluster, eh)) +
    expand_limits(y = c(0, 0.5)) +
    geom_violin(aes(fill = cluster)) +
    geom_boxplot(width = 0.2) +
    scale_fill_manual(
      "Cluster EH", values = colour_set[c(1, 2, 4)],
      na.translate = FALSE, drop = FALSE
    ) +
    remove_x_axis +
    labs(y = "Diversity") +
    guides(fill = FALSE)

  point_plots$cluster_eh <- ggplot() +
    scale_x_continuous(
      breaks = seq(0, max_pos_mb, by = 10), expand = c(0.01, 0.01)
    ) +
    expand_limits(x = c(0, max_pos_mb), y = c(0, 0.5)) +
    geom_point(
      aes(
        ehs_by_cluster[[chrom]]$pos_mb, ehs_by_cluster[[chrom]]$eh,
        colour = ehs_by_cluster[[chrom]]$cluster
      ), size = point_size
    ) +
    scale_colour_manual(
      "Cluster EH",
      values = colour_set[c(1, 2, 4)], na.translate = FALSE, drop = FALSE
    ) +
    plot_landmarks(0.5) +
    remove_x_axis +
    labs(y = "Diversity") +
    guides(linetype = FALSE)

  ## 3
  dist_plots$josts_d <- ggplot(
      snp_data[[chrom]], aes(josts_d_class, josts_d)
    ) +
    expand_limits(y = c(0, 1)) +
    geom_violin(aes(fill = josts_d_class)) +
    geom_boxplot(width = 0.2) +
    expand_limits(y = c(0, 1)) +
    scale_fill_manual(
      "Jost's D", values = colour_set[c(15, 17, 19)],
      na.translate = FALSE, drop = FALSE
    ) +
    remove_x_axis +
    labs(y = "Jost's D") +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    guides(fill = FALSE)

  point_plots$josts_d <- ggplot() +
    scale_x_continuous(
      breaks = seq(0, max_pos_mb, by = 10), expand = c(0.01, 0.01)
    ) +
    expand_limits(x = c(0, max_pos_mb), y = c(0, 1)) +
    geom_point(
      aes(
        snp_data[[chrom]]$pos_mb, snp_data[[chrom]]$josts_d,
        colour = snp_data[[chrom]]$josts_d_class,
      ), size = point_size
    ) +
    scale_colour_manual(
      "Jost's D Class", values = colour_set[c(15, 17, 19)],
      na.translate = FALSE, drop = FALSE
    ) +
    plot_landmarks(1) +
    remove_x_axis +
    labs(y = "Jost's D") +
    guides(linetype = FALSE)

  ## 4
  dist_plots$phys_pos_vs_gen_pos <- ggplot(
      snp_data[[chrom]], aes(chrom, pos_cm)
    ) +
    expand_limits(y = c(0, max_pos_cm)) +
    geom_violin(fill = colour_set[3]) +
    geom_boxplot(width = 0.2) +
    labs(x = "", y = "Position in cM")

  odi_present <- which(
    levels(snp_data[[chrom]]$odi) %in% as.character(snp_data[[chrom]]$odi)
  )

  point_plots$phys_pos_vs_gen_pos <- ggplot() +
    expand_limits(x = c(0, max_pos_mb), y = c(0, max_pos_cm)) +
    scale_x_continuous(
      breaks = seq(0, max_pos_mb, by = 10), expand = c(0.01, 0.01)
    ) +
    geom_point(
      aes(
        snp_data[[chrom]]$pos_mb , snp_data[[chrom]]$pos_cm,
        colour = snp_data[[chrom]]$odi
      ), size = point_size
    ) +
    scale_colour_manual(
      "Order\nDifference\nIntervals", values = colours_intervals[odi_present]
    ) +
    plot_landmarks(max_pos_cm) +
    labs(x = "Position in Mb", y = "Position in cM") +
    guides(linetype = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0))

  plot_types <- c("mjafs", "cluster_eh", "josts_d", "phys_pos_vs_gen_pos")

  lapply(
    plot_types, function (type) {
      list(
        dist_plots[[type]],
        point_plots[[type]] +
          theme(legend.position = "none") +
          remove_y_axis,
        as_ggplot(get_legend(point_plots[[type]]))
      )
    }
  ) %>% unlist(recursive = FALSE)
}, simplify = FALSE)

# ################################################################################
# # plot all plots

chrom_1A <- c(1:3, 13:15, 25:27)
group_1 <- lapply(seq(0, 9, 3), function (num) chrom_1A + num) %>% unlist()
order <- lapply(seq(0, 216, 36), function (num) group_1 + num) %>% unlist()

# plot the matrix
png(
  file.path("results", "phys_pos_vs_mjaf_cluster_mjaf_josts_d_and_gen_pos.png"),
  family = "Times New Roman", width = 4960, height = 7016, pointsize = 5
)
grid.arrange(
  grobs = (plots %>% unlist(recursive = FALSE))[order], nrow = 28, ncol = 9,
  widths = rep(c(6, 41, 3), 3)
)
dev.off()

# ################################################################################
# # plot just phys pos vs gen pos graphs

mod_plots <- lapply(chroms, function (chrom) {
  lapply(c(10, 11), function (index) {
    new_plot <- ggplot()
    if (chrom == "7A") {
      new_plot <- plots[[chrom]][[index]]

    } else if (grepl("A", chrom)) {
      new_plot <- plots[[chrom]][[index]] +
        remove_x_axis

    } else if (grepl("7", chrom)) {
      new_plot <- plots[[chrom]][[index]] +
        remove_y_axis
    } else {
      new_plot <- plots[[chrom]][[index]] +
        remove_y_axis +
        remove_x_axis
    }
    if (index == 10) {
      new_plot <- new_plot +
        labs(title = chrom)
    } else {
      new_plot <- new_plot +
        labs(title = chrom) +
        theme(plot.title = element_text(color = "white"))
    }
    new_plot
  })
})

# plot the matrix
png(
  file.path("results", "phys_pos_vs_gen_pos_v2.png"),
  family = "Times New Roman", width = 2480, height = 3508, pointsize = 5
)
grid.arrange(
  grobs = c(
    mod_plots %>% unlist(recursive = FALSE),
    (plots %>% unlist(recursive = FALSE))[24]
  ),
  layout_matrix = matrix(nrow = 6, 1:42) %>% rbind(rep(43, 7)) %>% t(),
  widths = c(rep(c(6, 44), 3), 5)
)
dev.off()

################################################################################
# plot zoomed windows

window_ranges <- read_csv(file.path(intermediate, "windows.csv"))

zoomed_plots <- lapply(chroms, function (chrom) {
  print(chrom)

  restore_y_axis <- theme(
    axis.title.y = element_text(angle = 90),
    axis.text.y = element_text(hjust = 1),
    axis.ticks.y = element_line(colour = "grey50")
  )

  restore_x_axis <- theme(
    axis.title.x = element_text(),
    axis.text.x = element_text(),
    axis.ticks.x = element_line(colour = "grey50")
  )

  base <- function (index, row) {
    plots[[chrom]][[index]] +
      xlim(
        chrom_window_ranges[row, ]$start, chrom_window_ranges[row, ]$end
      ) +
      expand_limits(y = c(0, 1)) +
      restore_y_axis
  }

  chrom_window_ranges <- window_ranges[which(window_ranges$chrom == chrom), ]

  if (nrow(chrom_window_ranges)) {
    lapply(1:nrow(chrom_window_ranges), function (row) {
      range <- str_c(
        chrom, "_",
        chrom_window_ranges[row, ]$start, "_",
        chrom_window_ranges[row, ]$end, "_",
        chrom_window_ranges[row, ]$genes
      )
      print(range)

      which_markers <- which(
        (snp_data[[chrom]]$pos_mb >= chrom_window_ranges[row, ]$start) &
        (snp_data[[chrom]]$pos_mb <= chrom_window_ranges[row, ]$end)
      )
      cm_range <- range(snp_data[[chrom]]$pos_cm[which_markers])

      zoomed_plots <- lapply(seq(2, 11, 3), function (index) {
        if (index == 2) {
          list(
            base(index, row) +
              scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)),
            plots[[chrom]][[index + 1]]
          )
        } else if (index == 5) {
          list(
            base(index, row),
            plots[[chrom]][[index + 1]]
          )
        } else if (index == 8) {
          list(
            base(index, row) +
              labs(x = "Position in Mb") +
              scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)),
            plots[[chrom]][[index + 1]]
          )
        } else {
          list(
            base(index, row) +
              restore_x_axis +
              scale_y_continuous(limits = cm_range),
            plots[[chrom]][[index + 1]]
          )
        }
      }) %>% unlist(recursive = FALSE)

      png(
        file.path(zoomed_marker_plots, str_c(range, ".png")),
        family = "Times New Roman", width = 270, height = 345,
        units = "mm", res = 96
      )
      grid.arrange(
        grobs = zoomed_plots, nrow = 4, ncol = 2, widths = c(42, 8)
      )
      dev.off()

      zoomed_plots
    }) %>% unlist(recursive = FALSE)
  }
})

# # zoomed plots just eh
# zoomed_plots <- lapply(chroms, function (chrom) {
#   print(chrom)

#   base <- function (index, row) {
#     plots[[chrom]][[index]] +
#       xlim(
#         chrom_window_ranges[row, ]$start, chrom_window_ranges[row, ]$end
#       ) +
#       expand_limits(y = c(0, 1)) +
#       theme(
#         axis.title.y = element_text(angle = 90),
#         axis.text.y = element_text(hjust = 1),
#         axis.ticks.y = element_line(colour = "grey50")
#       )
#   }

#   chrom_window_ranges <- window_ranges[which(window_ranges$chrom == chrom), ]

#   if (nrow(chrom_window_ranges)) {
#     lapply(1:nrow(chrom_window_ranges), function (row) {
#       range <- str_c(
#         chrom, "_",
#         chrom_window_ranges[row, ]$start, "_",
#         chrom_window_ranges[row, ]$end, "_",
#         chrom_window_ranges[row, ]$genes
#       )
#       print(range)

#       index <- 2
#       zoomed_plots <-  list(
#         base(index, row),
#         plots[[chrom]][[index + 1]]
#       )

#       png(
#         file.path(zoomed_marker_plots, str_c("just_mjaf_", range, ".png")),
#         family = "Times New Roman", width = 270, height = 90,
#         units = "mm", res = 96
#       )
#       grid.arrange(
#         grobs = zoomed_plots, nrow = 1, ncol = 2, widths = c(42, 8)
#       )
#       dev.off()

#       zoomed_plots
#     }) %>% unlist(recursive = FALSE)
#   }
# })

# # zoomed plots just maps
# zoomed_plots <- lapply(chroms, function (chrom) {
#   print(chrom)

#   base <- function (index, row) {
#     plots[[chrom]][[index]] +
#       xlim(
#         chrom_window_ranges[row, ]$start, chrom_window_ranges[row, ]$end
#       ) +
#       expand_limits(y = c(0, 1)) +
#       theme(
#         axis.title.y = element_text(angle = 90),
#         axis.text.y = element_text(hjust = 1),
#         axis.ticks.y = element_line(colour = "grey50")
#       )
#   }

#   chrom_window_ranges <- window_ranges[which(window_ranges$chrom == chrom), ]

#   if (nrow(chrom_window_ranges)) {
#     lapply(1:nrow(chrom_window_ranges), function (row) {
#       range <- str_c(
#         chrom, "_",
#         chrom_window_ranges[row, ]$start, "_",
#         chrom_window_ranges[row, ]$end, "_",
#         chrom_window_ranges[row, ]$genes
#       )
#       print(range)

#       which_markers <- which(
#         (snp_data[[chrom]]$pos_mb >= chrom_window_ranges[row, ]$start) &
#         (snp_data[[chrom]]$pos_mb <= chrom_window_ranges[row, ]$end)
#       )
#       cm_range <- range(snp_data[[chrom]]$pos_cm[which_markers])

#       index <- 11
#       zoomed_plots <-  list(
#         base(index, row) +
#           theme(
#             axis.title.x = element_text(),
#             axis.text.x = element_text(),
#             axis.ticks.x = element_line(colour = "grey50")
#           ) +
#           scale_y_continuous(limits = cm_range),
#         plots[[chrom]][[index + 1]]
#       )

#       png(
#         file.path(zoomed_marker_plots, str_c("just_maps_", range, ".png")),
#         family = "Times New Roman", width = 270, height = 90,
#         units = "mm", res = 96
#       )
#       grid.arrange(
#         grobs = zoomed_plots, nrow = 1, ncol = 2, widths = c(42, 8)
#       )
#       dev.off()

#       zoomed_plots
#     }) %>% unlist(recursive = FALSE)
#   }
# })