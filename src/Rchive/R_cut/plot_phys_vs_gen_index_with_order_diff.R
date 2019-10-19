source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(pgda, "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_point", "labs",  "scale_colour_manual", 
  "theme", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(stringr, "str_c", "str_wrap")
import::from(tibble, "tibble")

################################################################################
# import data

phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

################################################################################
# caclulate each markers order difference between gen and phys maps and assign
# it to an interval

gen_to_phys_order <- match(phys_data$snp$id, gen_data$snp$id)
order_diffs <- (gen_to_phys_order - 1:length(gen_to_phys_order)) %>% abs()

order_diff_quantiles <- c(
  "0%" = -1,
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
# collect the data

snp_data <- tibble(
  phys_id = phys_data$snp$id,
  gen_id = gen_data$snp$id,
  chrom = phys_data$snp$chrom,
  phys_pos_mb = phys_data$snp$pos / 1e6,
  gen_pos_cm = gen_data$snp$pos[match(phys_data$snp$id, gen_data$snp$id)] / 100,
  odi = order_diff_intervals
)

################################################################################
# calc the lengths of the different genomes and homoeologous sets

max_markers <- by(snp_data, snp_data$chrom, nrow) %>% max_lengths()

################################################################################
# plot the data

plots <- by(snp_data, snp_data$chrom, function (chrom_data) {
  chrom <- chrom_data$chrom[1]

  gen_to_phys_order <- match(chrom_data$phys_id, chrom_data$gen_id)

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
        1:length(gen_to_phys_order), gen_to_phys_order,
        colour = chrom_data$odi
      ), size = 0.3
    ) +
    # labs(colour = levels) +
    scale_colour_manual(name = "Order Difference", values = colours_intervals)
})

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  xlab = "Pseudo-Chromosomal Index", ylab = "Genetic Index",
  title = str_wrap(
    str_c(
      "Pseudo-Chromosomal vs Genetic Order with Markers Coloured by ",
      "Magnitude of Difference in Order Between Maps"
    ), 
    width = 70
  ),
  legend = c(2, 1)
)

# plot the matrix
png(
  file.path("results", "phys_vs_gen_index_with_order_diff.png"),
  family = "Times New Roman", width = 360, height = 640, pointsize = 5,
  units = "mm", res = 96
)
plots_matrix + theme(legend.position = "bottom")
dev.off()

# summary(order_diffs)
# quantile(order_diffs, c(0.5, 0.75, 0.875, 0.9375, 0.96875, 0.984375, 0.992187))

