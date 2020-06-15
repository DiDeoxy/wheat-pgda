source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(tidyr, "gather")
import::from(pgda, "calc_eh", "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_smooth", "labs",  "scale_colour_gradientn", 
  "scale_size_continuous", "theme", "unit", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(readr, "read_rds")
import::from(Rfast, "rowMaxs")
import::from(scrime, "rowTables")
import::from(stringr, "str_c")
import::from(tibble, "as_tibble", "tibble")

# load the data from the gds object
wheat_data <- snpgds_parse(phys_gds)

# calc the lengths of the different genomes and homoeologous sets
max_lengths <- span_by_chrom(
  wheat_data$snp$chrom, wheat_data$snp$pos
) %>% max_lengths() / 1e6

# recode the genotypes
wheat_data$genotypes <- replace(wheat_data$genotypes, wheat_data$geno == 3, NA)

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

# calc each markers mjaf by each cluster group
snp_data <- tibble(
  chrom = wheat_data$snp$chrom,
  pos_mb = wheat_data$snp$pos / 1e6,
  josts_d = read_rds(josts_d)
)

coding <- c(0, 2)
mja <- rowMaxs(rowTables(wheat_data$genotypes, coding))
mjafs_by_pop <- lapply(genos, function (geno) {
  geno_counts <- rowTables(geno, coding)
  max_genos <- geno_counts[cbind(seq_along(mja), mja)]
  max_genos / rowSums(geno_counts)
}) %>% do.call(cbind, .)

snp_data <- cbind(snp_data, mjafs_by_pop) %>% as_tibble()

top_josts_d <- which(snp_data$josts_d >= quantile(snp_data$josts_d, 0.75))
snp_data <- snp_data[top_josts_d, ]

snp_data_long <- snp_data %>% gather(group, mjaf, -c(chrom, pos_mb, josts_d))

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
      ylim(0, 1) +
      geom_point(
        aes(pos_mb, mjaf, colour = group), size = 0.5, alpha = 0.5
      ) +
      scale_colour_manual(
        "Cluster Group", values = colours_comps_genes,
        labels = c("CHRS", "CHRW", "CSWS"),
      )
  }
)

empty_plot <- ggplot() +
  xlim(0, max_lengths[["D"]]) + ylim(0, 1) + geom_point(aes(0.5, 10), alpha = 0)

plots <- append(plots, list("6D" = empty_plot), 17)

plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3,
  xlab = "Position in Mb",
  ylab = "Mjaor Allele Frequency",
  xAxisLabels = c("A", "B", "D"),
  yAxisLabels = 1:7,
  legend = c(1, 1)
)

# # plot the matrix
png(
  file.path("results", "phys_pos_vs_mjaf_by_group.png"),
  family = "Times New Roman", width = 165, height = 208, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()