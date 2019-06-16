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

phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

snp_data <- tibble(
  phys_id = phys_data$snp$id,
  gen_id = gen_data$snp$id,
  chrom = phys_data$snp$chrom,
  phys_pos_mb = phys_data$snp$pos / 1e6,
  gen_pos_cm = gen_data$snp$pos / 100
)

# calc the lengths of the different genomes and homoeologous sets
max_markers <- by(snp_data, snp_data$chrom, nrow) %>% max_lengths()

# create a function for making a gradient of colours
levels <- c(
  "0-7", "8-26", "27-72", "73-126", "127-173", "174-298", "299-436", "437-1248"
)
names(colours_order_diff) <- levels

# use to hold all diffs for calcing stats from
order_diffs <- c()

plots <- by(snp_data, snp_data$chrom, function (chrom_data) {
  chrom <- chrom_data$chrom[1]

  gen_to_phys_order <- match(chrom_data$phys_id, chrom_data$gen_id)
  order_diff <- (
    gen_to_phys_order - 1:length(gen_to_phys_order)
  ) %>% abs()

  order_diffs <<- c(order_diffs, order_diff)

  order_diff_intervals <- cut(
    order_diff, c(-1, 7, 26, 72, 126, 173, 298, 436, 1248), levels
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
        1:length(gen_to_phys_order), gen_to_phys_order,
        colour = order_diff_intervals
      ), size = 0.3
    ) +
    # labs(colour = levels) +
    scale_colour_manual(name = "Order Difference", values = colours_order_diff)
})

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  xlab = "Pseudo-Chromosomal Index", ylab = "Genetic Index",
  # title = str_wrap(
  #   str_c(
  #     "Pseudo-Chromosomal vs Genetic Order with Markers Coloured by ",
  #     "Magnitude of Difference in Order Between Maps"
  #   ), 
  #   width = 70
  # ),
  legend = c(2, 1)
)

# plot the matrix
png(
  file.path("results", "phys_vs_gen_index_with_order_diff.png"),
  family = "Times New Roman", width = 165, height = 208, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom")
dev.off()

summary(order_diffs)
quantile(order_diffs, c(0.5, 0.75, 0.875, 0.9375, 0.96875, 0.984375, 0.992187))

