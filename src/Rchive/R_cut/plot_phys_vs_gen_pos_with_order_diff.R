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
max_phys_lengths <- span_by_chrom(
  phys_data$snp$chrom, phys_data$snp$pos
) %>% max_lengths() / 1e6
max_gen_lengths <- span_by_chrom(
  gen_data$snp$chrom, gen_data$snp$pos
) %>% max_lengths() / 100

# create a function for making a gradient of colours
levels <- c(
  "0-7", "8-26", "27-72", "73-126", "127-173", "174-298", "299-436", "437-1248"
)
names(colours_order_diff) <- levels

plots <- by(snp_data, snp_data$chrom, function (chrom_data) {
  chrom <- chrom_data$chrom[1]

  gen_to_phys_order <- match(chrom_data$phys_id, chrom_data$gen_id)
  order_diff <- (
    gen_to_phys_order - 1:length(gen_to_phys_order)
  ) %>% abs()

  order_diff_intervals <- cut(
    order_diff, c(-1, 7, 26, 72, 126, 173, 298, 436, 1248), levels
  )

  ggplot() +
    xlim(
      0,
      max_phys_lengths[[
        ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
      ]]
    ) +
    ylim(
      0,
      max_gen_lengths[[
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
        chrom_data$phys_pos_mb, chrom_data$gen_pos_cm[gen_to_phys_order],
        colour = order_diff_intervals
      ), size = 1.5
    ) +
    labs(colour = levels) +
    scale_colour_manual(name = "Order Difference", values = colours_order_diff)
})

# plots_matrix <- ggmatrix(
#   plots, nrow = 7, ncol = 3, xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
#   xlab = "Pseudo-Chromosomal Position in Mb", ylab = "Genetic Position in cM",
#   # title = str_wrap(
#   #   str_c(
#   #     "Pseudo-Chromosomal vs Genetic Position with Markers Coloured by ",
#   #     "Magnitude of Difference in Order Between Maps"
#   #   ), width = 70
#   # ),
#   legend = c(2, 1)
# )

plots_matrix_pos <- ggmatrix(
  plots[c(6, 10)], nrow = 2, ncol = 1, yAxisLabels = c("2D", "4A"),
  xlab = "Position in Mb", ylab = "Position in cM",
  legend = c(2, 1),
  title = str_c(
    "Pseudo-chromosomal vs Genetic position with markers\n",
    "coloured by magnitude of order difference"
  )
)

# plot the matrix
png(
  file.path("results", "2D_4A_phys_vs_gen_pos_with_order_diff.png"),
  family = "Times New Roman", width = 160, height = 192, pointsize = 5,
  units = "mm", res = 300)
plots_matrix_pos + theme(legend.position = "bottom")
dev.off()