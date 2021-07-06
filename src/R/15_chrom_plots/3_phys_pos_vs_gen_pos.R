source("wheat-pgda/src/R/15_chrom_plots/1_base_data.R")
library(ggplot2)
import::from(ggpubr, "as_ggplot", "get_legend")
import::from(ggrepel, "geom_label_repel")
import::from(gridExtra, "grid.arrange")

line_size <- 1
point_size <- 3
text_size <- 9
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

  dist_plot <- ggplot(
      snp_data[[chrom]], aes(chrom, pos_cm)
    ) +
    expand_limits(y = c(0, max_pos_cm)) +
    geom_violin(fill = colour_set[3]) +
    geom_boxplot(width = 0.2) +
    labs(x = "", y = "Position in cM", title = chrom) +
    theme(
      axis.text.x = element_text(size = text_size * 3, colour = "white"),
      axis.text.y = element_text(size = text_size * 3, angle = 90),
      # axis.title = element_blank(),
      plot.title = element_text(size = text_size * 3)
    )

  odi_present <- which(
    levels(snp_data[[chrom]]$odi) %in% as.character(snp_data[[chrom]]$odi)
  )

  point_plot <- ggplot() +
    expand_limits(x = c(0, max_pos_mb), y = c(0, max_pos_cm)) +
    scale_x_continuous(
      breaks = seq(0, max_pos_mb, by = 100), expand = c(0.01, 0.01)
    ) +
    geom_point(
      aes(
        snp_data[[chrom]]$pos_mb , snp_data[[chrom]]$pos_cm,
        colour = snp_data[[chrom]]$odi
      ), size = point_size
    ) +
    scale_colour_manual(
      "Order\nDifference", values = colours_intervals[odi_present]
    ) +
    plot_landmarks(max_pos_cm) +
    labs(x = "Position in Mb", y = "Position in cM", title = chrom) +
    guides(linetype = FALSE) +
    theme(
      axis.text.x = element_text(size = text_size * 3),
      # axis.title = element_blank(),
      plot.title = element_text(color = "white", size = text_size * 3),
      legend.position = "bottom",
      legend.text = element_text(size = text_size * 3),
      legend.title = element_text(size = text_size * 2.5)
    )

  if (chrom == "7A") {
    point_plot <- point_plot + remove_y_axis
  } else if (grepl("A", chrom)) {
    dist_plot <- dist_plot + remove_x_axis
    point_plot <- point_plot + remove_x_axis + remove_y_axis
  } else if (grepl("7", chrom)) {
    dist_plot <- dist_plot + remove_y_axis
    point_plot <- point_plot + remove_y_axis
  } else {
    dist_plot <- dist_plot + remove_x_axis + remove_y_axis
    point_plot <- point_plot + remove_x_axis + remove_y_axis
  }

  list(
    dist_plot,
    point_plot + theme(legend.position = "none"),
    as_ggplot(get_legend(point_plot))
  )

}, simplify = FALSE)

chrom_1A <- 1:2
group_1 <- lapply(seq(0, 6, 3), function (num) chrom_1A + num) %>% unlist()
order <- lapply(seq(0, 54, 9), function (num) group_1 + num) %>% unlist()

# plot the matrix
png(
  file.path("results", "phys_pos_vs_gen_pos.png"),
  family = "Times New Roman", width = 2480, height = 3508
)
grid.arrange(
  grobs = c(
    (plots %>% unlist(recursive = FALSE))[order],
    (plots %>% unlist(recursive = FALSE))[6]
  ),
  layout_matrix = matrix(nrow = 6, 1:42) %>% cbind(rep(43, 6)) %>% t(),
  widths = rep(c(6, 44), 3),
  heights = c(rep(10, 7), 2)
)
dev.off()