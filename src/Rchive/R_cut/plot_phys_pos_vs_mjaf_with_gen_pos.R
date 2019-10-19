source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(pgda, "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_point", "labs",  "scale_colour_gradientn", 
  "scale_size_continuous", "theme", "unit", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(Rfast, "rowMaxs")
import::from(scrime, "rowTables")
import::from(stringr, "str_c")
import::from(tibble, "tibble")

# load the data from the gds object
phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

# recode the genotypes and count
genos <- replace(phys_data$genotypes, phys_data$geno == 3, NA)
allele_counts <- rowTables(genos, c(0, 2))

# make a tibble with the relevant data
snp_data <- tibble(
  chrom = phys_data$snp$chrom,
  phys_pos_mb = phys_data$snp$pos / 1e6,
  gen_pos_cm = gen_data$snp$pos[match(phys_data$snp$id, gen_data$snp$id)] / 100,
  mjafs = rowMaxs(allele_counts,  value = TRUE) / rowSums(allele_counts)
)

# create a function for making a gradient of colours
# colour_gradient <- colorRampPalette(colour_set[c(1, 5, 3, 2, 4)])(100)
colour_gradient <- gray.colors(100, start = 0.8, end = 0.1)

# calc the lengths of the different genomes and homoeologous sets
max_phys_lengths <- span_by_chrom(
  phys_data$snp$chrom, phys_data$snp$pos
) %>% max_lengths() / 1e6
max_gen_lengths <- span_by_chrom(
  gen_data$snp$chrom, gen_data$snp$pos
) %>% max_lengths() / 100

# create plots of phys position vs gen pos
plots <- by(snp_data, snp_data$chrom,
  function (chrom_data) {
    chrom <- chrom_data$chrom[1]
    chrom_data %>%
      ggplot() +
      xlim(
        0,
        max_phys_lengths[[
          ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
        ]]
      ) +
      # ylim(
      #   0,
      #   max_gen_lengths[[
      #     ifelse(grepl("1", chrom), "one", 
      #       ifelse(grepl("2", chrom), "two",
      #         ifelse(grepl("3", chrom), "three",
      #           ifelse(grepl("4", chrom), "four",
      #             ifelse(grepl("5", chrom), "five",
      #               ifelse(grepl("6", chrom), "six", "seven")
      #             )
      #           )
      #         )
      #       )
      #     )
      #   ]]
      # ) +
      ylim(0.5, 1) +
      geom_point(
        aes(
          phys_pos_mb, mjafs, colour = (gen_pos_cm / max(max_gen_lengths) * 100)
        ), size = 0.5
      ) +
      scale_colour_gradientn(
         name = "Percent Max cM Position", colours = colour_gradient
      ) +
      theme(legend.key.size = unit(15, "points"))
      # geom_point(
      #   aes(haplo_phys, haplo_gen),
      #   colour = "black",
      #   shape = 1, size = 0.75
      # )
  }
)

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3,
  xlab = "Position in Mb", ylab = "Major Allele Frequency",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_c(
    "Marker Pseudo-Chromosomal vs Major Allele Frequency by\n",
    "Chromosome Coloured by Markers' Percent of Global\n",
    "Maximum Genetic Position"
  ),
  legend = c(5, 1)
)

# plot the matrix
png(
  file.path("results", "phys_vs_mjaf_with_gen_pos.png"),
  family = "Times New Roman", width = 160, height = 192, pointsize = 5,
  units = "mm", res = 192)
plots_matrix + theme(legend.position = "bottom")
dev.off()