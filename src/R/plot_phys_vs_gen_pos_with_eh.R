# load the data from the gds object
wheat_data_phys <- parse_gds(phys_gds)
wheat_data_gen <- parse_gds(gen_gds)

snp_phys_order <- match(wheat_data_phys$snp$id, wheat_data_gen$snp$id)

# make a tibble with the relevant data
phys_gen_snp_pos <- tibble(
  chrom = wheat_data_phys$snp$chrom, phys = wheat_data_phys$snp$pos_mb,
  gen = wheat_data_gen$snp$pos[snp_phys_order] / 100,
  eh = calc_eh(wheat_data_phys$genotypes)
)

# # identify those snps within the extended haplotypes
# snp_index_1A <- which(wheat_data_phys$snp$chrom == "1A")
# haplo_id_1A <- wheat_data_phys$snp$id[snp_index_1A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_1A] > 70 &
#     wheat_data_phys$snp$pos_mb[snp_index_1A] < 300)
# ]
# snp_index_2A <- which(wheat_data_phys$snp$chrom == "2A")
# haplo_id_2A <- wheat_data_phys$snp$id[snp_index_2A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_2A] > 210 &
#     wheat_data_phys$snp$pos_mb[snp_index_2A] < 470)
# ]
# snp_index_4A <- which(wheat_data_phys$snp$chrom == "4A")
# haplo_id_4A <- wheat_data_phys$snp$id[snp_index_4A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_4A] > 230 &
#     wheat_data_phys$snp$pos_mb[snp_index_4A] < 460)
# ]
# snp_index_5B <- which(wheat_data_phys$snp$chrom == "5B")
# haplo_id_5B <- wheat_data_phys$snp$id[snp_index_5B][
#   which(wheat_data_phys$snp$pos_mb[snp_index_5B] > 110 &
#     wheat_data_phys$snp$pos_mb[snp_index_5B] < 210)
# ]
# snp_index_6A <- which(wheat_data_phys$snp$chrom == "6A")
# haplo_id_6A <- wheat_data_phys$snp$id[snp_index_6A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_6A] > 170 &
#     wheat_data_phys$snp$pos_mb[snp_index_6A] < 445)
# ]
# snp_index_6B <- which(wheat_data_phys$snp$chrom == "6B")
# haplo_id_6B <- wheat_data_phys$snp$id[snp_index_6B][
#   which(wheat_data_phys$snp$pos_mb[snp_index_6B] > 250 &
#     wheat_data_phys$snp$pos_mb[snp_index_6B] < 380)
# ]
# snp_index_7A <- which(wheat_data_phys$snp$chrom == "7A")
# haplo_id_7A <- wheat_data_phys$snp$id[snp_index_7A][
#   which(wheat_data_phys$snp$pos_mb[snp_index_7A] > 310 &
#     wheat_data_phys$snp$pos_mb[snp_index_7A] < 445)
# ]
# haplo_ids <- c(
#   haplo_id_1A, haplo_id_2A, haplo_id_4A, haplo_id_5B, haplo_id_6A, haplo_id_6B,
#   haplo_id_7A
# )
# haplo_index_snps_phys <- match(haplo_ids, wheat_data_phys$snp$id)
# haplo_index_snps_gen <- match(haplo_ids, wheat_data_gen$snp$id)

# # create a column in the snp_data set that contains D values for only
# # those markers in the extended haplotypes, NA for all else
# phys_gen_snp_pos <- phys_gen_snp_pos %>% mutate(haplo_phys = phys)
# phys_gen_snp_pos$haplo_phys[-haplo_index_snps_phys] <- NA
# phys_gen_snp_pos <- phys_gen_snp_pos %>% mutate(haplo_gen = gen)
# phys_gen_snp_pos$haplo_gen[-haplo_index_snps_gen] <- NA

# allows application of same colour to each set of chromosomes
chroms_order <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>% 
  t() %>% as.vector()
colour_order <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3),
rep(6, 3), rep(7, 3))

# create a function for making a gradient of colours
colour_gradient <- colorRampPalette(c("Red", "Green", "Blue"))

# create plots of phys position vs gen pos
plots <- by(phys_gen_snp_pos, phys_gen_snp_pos$chrom,
  function (data_chrom) {
    chrom <- data_chrom$chrom[1]
    data_chrom %>%
      ggplot() +
      xlim(
        0,
        wheat_data_phys$max_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ]
      ) +
      ylim(
        0,
        wheat_data_gen$max_lengths[
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
        ]
      ) +
      geom_point(aes(phys, gen, colour = eh), size = 0.5) +
      labs(
        colour = "Expected Heterozygosity", size = "Linkage Disequilibrium"
      ) +
      scale_colour_gradientn(colours = colour_gradient(100)) +
      scale_size_continuous() +
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
    "Comparison of Order by Position between Location and Genetic Maps"
  ),
  legend = c(1, 1)
)

# plot the matrix
png(
  file.path("results", "chrom_marker_phys_vs_gen_pos_and_eh.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300)
plots_matrix + theme(legend.position = "bottom")
dev.off()