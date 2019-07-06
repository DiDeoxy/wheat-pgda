source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(pgda, "calc_eh", "max_lengths", "snpgds_parse", "span_by_chrom")
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

# # identify those snps within the extended haplotypes
# snp_index_1A <- which(phys_data$snp$chrom == "1A")
# haplo_id_1A <- phys_data$snp$id[snp_index_1A][
#   which(phys_data$snp$pos_mb[snp_index_1A] > 70 &
#     phys_data$snp$pos_mb[snp_index_1A] < 300)
# ]
# snp_index_2A <- which(phys_data$snp$chrom == "2A")
# haplo_id_2A <- phys_data$snp$id[snp_index_2A][
#   which(phys_data$snp$pos_mb[snp_index_2A] > 210 &
#     phys_data$snp$pos_mb[snp_index_2A] < 470)
# ]
# snp_index_4A <- which(phys_data$snp$chrom == "4A")
# haplo_id_4A <- phys_data$snp$id[snp_index_4A][
#   which(phys_data$snp$pos_mb[snp_index_4A] > 230 &
#     phys_data$snp$pos_mb[snp_index_4A] < 460)
# ]
# snp_index_5B <- which(phys_data$snp$chrom == "5B")
# haplo_id_5B <- phys_data$snp$id[snp_index_5B][
#   which(phys_data$snp$pos_mb[snp_index_5B] > 110 &
#     phys_data$snp$pos_mb[snp_index_5B] < 210)
# ]
# snp_index_6A <- which(phys_data$snp$chrom == "6A")
# haplo_id_6A <- phys_data$snp$id[snp_index_6A][
#   which(phys_data$snp$pos_mb[snp_index_6A] > 170 &
#     phys_data$snp$pos_mb[snp_index_6A] < 445)
# ]
# snp_index_6B <- which(phys_data$snp$chrom == "6B")
# haplo_id_6B <- phys_data$snp$id[snp_index_6B][
#   which(phys_data$snp$pos_mb[snp_index_6B] > 250 &
#     phys_data$snp$pos_mb[snp_index_6B] < 380)
# ]
# snp_index_7A <- which(phys_data$snp$chrom == "7A")
# haplo_id_7A <- phys_data$snp$id[snp_index_7A][
#   which(phys_data$snp$pos_mb[snp_index_7A] > 310 &
#     phys_data$snp$pos_mb[snp_index_7A] < 445)
# ]
# haplo_ids <- c(
#   haplo_id_1A, haplo_id_2A, haplo_id_4A, haplo_id_5B, haplo_id_6A, haplo_id_6B,
#   haplo_id_7A
# )
# haplo_index_snps_phys <- match(haplo_ids, phys_data$snp$id)
# haplo_index_snps_gen <- match(haplo_ids, gen_data$snp$id)

# # create a column in the snp_data set that contains D values for only
# # those markers in the extended haplotypes, NA for all else
# snp_data <- snp_data %>% mutate(haplo_phys = phys)
# snp_data$haplo_phys[-haplo_index_snps_phys] <- NA
# snp_data <- snp_data %>% mutate(haplo_gen = gen)
# snp_data$haplo_gen[-haplo_index_snps_gen] <- NA

# allows application of same colour to each set of chromosomes
chroms_order <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()
colour_order <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3),
rep(6, 3), rep(7, 3))

# create a function for making a gradient of colours
colour_gradient <- colorRampPalette(colour_set[c(1, 5, 3, 2, 4)])(100)

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
      geom_point(aes(phys_pos_mb, mjafs, colour = gen_pos_cm), size = 0.5) +
      scale_colour_gradientn(
         name = "Major Allele Frequency", colours = colour_gradient
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
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Position in cM",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_c(
    "Marker Pseudo-Chromosomal vs Genetic Position by\n",
    "Chromosome Coloured by Major Allele Frequency"
  ),
  legend = c(1, 1)
)

# plot the matrix
png(
  file.path("results", "phys_vs_gen_pos_with_mjaf.png"),
  family = "Times New Roman", width = 200, height = 240, pointsize = 5,
  units = "mm", res = 192)
plots_matrix + theme(legend.position = "bottom")
dev.off()