source("wheat-pgda/src/R/15_chrom_plots/1_base_data.R")
library(ggplot2)
import::from(ggpubr, "as_ggplot", "get_legend")
import::from(ggrepel, "geom_label_repel")
import::from(gridExtra, "grid.arrange")
import::from(magrittr, "%>%")

line_size <- 1
point_size <- 2
text_size <- 5
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