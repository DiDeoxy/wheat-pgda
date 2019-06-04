source(file.path("src", "R", "file_paths.R"))
library(tidyverse)
library(GGally)
library(pgda)
library(parallel)
library(mgcv)
library(RColorBrewer)

phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

snp_data <- tibble(
  phys_id = phys_data$snp$id, gen_id = gen_data$snp$id,
  chrom = phys_data$snp$chrom,
  phys_pos_mb = phys_data$snp$pos / 1e6,
  gen_pos_cm = gen_data$snp$pos / 100
)

# create a function for making a gradient of colours
levels <- c(
  "0-10", "11-25", "26-50", "51-100", "101-250", "251-500", "501-1000",
  "1001-2500"
)
colour_levels <- brewer.pal(n = 8, name = "Set1")[c(2, 1, 3:8)]
names(colour_levels) <- levels

order_diffs <- c()

# calc the lengths of the different genomes and homoeologous sets
max_markers <- by(snp_data, snp_data$chrom, nrow) %>% max_lengths()

# create a function for making a gradient of colours
levels <- c(
  "0-4", "5-10", "11-15", "16-20", "21-50", "51-100", "101-500"
)
colour_levels <- brewer.pal(n = 8, name = "Set1")[c(2, 3, 6, 5, 1, 4, 7)]
names(colour_levels) <- levels

plots <- by(snp_data, snp_data$chrom, function (chrom_data) {
  chrom <- chrom_data$chrom[1]

  relative_gen_order_forward <- match(chrom_data$phys_id, chrom_data$gen_id)
  relative_gen_order_reverse <- match(
    rev(chrom_data$phys_id), rev(chrom_data$gen_id)
  )

  relative_gen_order_forward_diff <- (
    relative_gen_order_forward - 1:length(relative_gen_order_forward)
  ) %>% abs()
  relative_gen_order_reverse_diff <- (
    relative_gen_order_reverse - 1:length(relative_gen_order_reverse)
  ) %>% abs() %>% rev()

  # order_diff <- (
  #   relative_gen_order_forward - 1:length(relative_gen_order_forward)
  # ) %>% abs()
  order_diff <- pmin(
    relative_gen_order_forward_diff, relative_gen_order_reverse_diff
  )

  order_diffs <<- c(order_diffs, order_diff)

  order_diff_intervals <- cut(
    order_diff, c(-1, 4, 10, 15, 20, 50, 100, 500), levels
  )

  ggplot() +
    xlim(
      0,
      max_markers[[
        ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
      ]]
    ) +
    ylim(
      0,
      max_markers[[
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
    ) +
    geom_point(
      aes(
        1:length(relative_gen_order_forward), relative_gen_order_forward,
        colour = order_diff_intervals
      ), size = 0.5
    ) +
    labs(colour = levels) +
    scale_colour_manual(name = "Min. Order Difference", values = colour_levels)
})

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Pseudo-Chromosomal Order",
  ylab = "Genetic Order", xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_wrap(
    str_c(
      "Pseudo-Chromosomal vs Genetic Order with Markers Coloured by ",
      "Magnitude of Minimum Difference in Order Between Maps"
    ), width = 50
  ),
  legend = c(2, 1)
)

# # plot the matrix
# png(
#   file.path("results", "phys_vs_gen_min_order_diff.png"),
#   family = "Times New Roman", width = 120, height = 267, pointsize = 5,
#   units = "mm", res = 300)
plots_matrix + theme(legend.position = "bottom")
# dev.off()

summary(order_diffs)