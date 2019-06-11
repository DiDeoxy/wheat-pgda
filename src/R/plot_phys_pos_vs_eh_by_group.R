source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(pgda, "calc_eh", "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_smooth", "labs",  "scale_colour_gradientn", 
  "scale_size_continuous", "theme", "unit", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(stringr, "str_c")
import::from(tibble, "tibble")

# load the data from the gds object
wheat_data <- snpgds_parse(phys_gds)

# calc the lengths of the different genomes and homoeologous sets
max_lengths <- span_by_chrom(
  wheat_data$snp$chrom, wheat_data$snp$pos
) %>% max_lengths() / 1e6

# get the cluster groups
cluster <- read_rds(hdbscan)$cluster
genos <- list(
  chrs = wheat_data$genotype[,
    which(wheat_data$sample$annot$pheno == "HRS" & cluster == 5)
  ],
  chrw = wheat_data$genotype[,
    which(wheat_data$sample$annot$pheno == "HRW" & cluster == 1)
  ],
  csws = wheat_data$genotype[,
    which(wheat_data$sample$annot$pheno == "SWS" & cluster == 2)
  ]
)

# make a tibble with the relevant data
snp_data <- tibble(
  chrom = wheat_data$snp$chrom,
  pos_mb = wheat_data$snp$pos / 1e6,
  CHRS = calc_eh(genos$chrs),
  CHRW = calc_eh(genos$chrw),
  CSWS = calc_eh(genos$csws)
)

snp_data_long <- snp_data %>% gather(group, eh, -c(chrom, pos_mb))

# create plots of phys position vs gen pos
plots <- by(snp_data_long, snp_data_long$chrom,
  function (chrom_data) {
    chrom <- chrom_data$chrom[1]
    chrom_data %>%
      ggplot() +
      xlim(
        0,
        max_lengths[[
          ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
        ]]
      ) +
      ylim(0, 0.5) +
      geom_smooth(
        aes(pos_mb, eh, colour = group), size = 0.75, se = FALSE,
        method = "loess", span = 0.1
      ) +
      geom_point(
        aes(pos_mb, eh, colour = group), size = 0.5, alpha = 0.1
      ) +
      scale_colour_manual(
        "Cluster Group", values = colours_comps_genes
      )
  }
)

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = "Expected Heterozygosity",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7, legend = c(1, 1)
)

# # plot the matrix
png(
  file.path("results", "phys_pos_vs_eh_by_group.png"),
  family = "Times New Roman", width = 165, height = 208, pointsize = 5,
  units = "mm", res = 300)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()