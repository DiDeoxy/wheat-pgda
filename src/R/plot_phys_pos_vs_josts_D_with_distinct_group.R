# load file paths, colours, and functions
source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "arrange", "mutate", "rowwise", "select")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2,
  "aes", "element_text", "ggplot", "geom_line", "geom_point", "geom_smooth",
  "guide_legend", "labs", "scale_colour_manual", "scale_x_continuous", "theme",
  "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(
  pgda, "calc_eh", "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom"
)
import::from(plyr, "rbind.fill")
import::from(readr, "read_rds", "write_csv")
import::from(Rfast, "rowMaxs")
import::from(scrime, "rowTables")
import::from(stringr, "str_c", "str_wrap")
import::from(tibble, "add_column", "add_row", "as_tibble", "tibble")
import::from(tidyr, "gather")

# load the data from the gds object
wheat_data <- snpgds_parse(phys_gds)

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
coding <- c(0, 2)
mja <- rowMaxs(rowTables(wheat_data$genotypes, coding))
mjafs_by_pop <- lapply(genos, function (geno) {
  geno_counts <- rowTables(geno, coding)
  max_genos <- geno_counts[cbind(seq_along(mja), mja)]
  max_genos / rowSums(geno_counts)
}) %>% do.call(cbind, .)

# add the genes positons to the regions table
snp_data <- wheat_data$snp %>%
  add_column(josts_d := read_rds(josts_d)) %>%
  cbind(mjafs_by_pop %>% round(4))
top_quartile <- snp_data$josts_d %>% quantile(0.75, na.rm = T)
snp_data <- snp_data %>%
  rowwise() %>%
  mutate(class =
    ifelse(
      josts_d >= top_quartile,
      c("CHRSD", "CHRWD", "CSWSD")[
        which.max(
          c(
            sum(abs(chrs - chrw), abs(chrs - csws)),
            sum(abs(chrw - chrs), abs(chrw - csws)),
            sum(abs(csws - chrs), abs(csws - chrw))
          )
        )
      ],
      "Indistinct"
    )
  )

# load the gene positions
pheno_genes <- load_genes(
  file.path(blast, "selected_pheno.csv"), base = 1
) %>% mutate(gene_type = "Phenotype Genes")
resi_genes <- load_genes(
  file.path(blast, "selected_resi.csv"), base = 1
) %>% mutate(gene_type = "Resistance Genes")

snp_data <- snp_data %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  mutate(pos_mb = pos / 1e6) %>%
  arrange(chrom, pos_mb) %>%
  mutate(type = pmin(gene_type, class, na.rm = TRUE)) %>%
  select(-c(gene_type, class, pos))

# calc the lengths of the different genomes and homoeologous sets
max_lengths <- span_by_chrom(
  wheat_data$snp$chrom, wheat_data$snp$pos
) %>% max_lengths() / 1e6

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
legend_title <- "Marker Types & Genes"
plots <- by(
  snp_data, snp_data$chrom, function(chrom_groups) {
    chrom <- chrom_groups$chrom[1]
    chrom_groups %>%
      ggplot() +
      ylim(0, 1) +
      xlim(
        0,
        max_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ]
      ) +
      geom_point(
        aes(pos_mb, josts_d, colour = type), shape = 16, size = 1
      ) +
      geom_smooth(
        aes(pos_mb, josts_d), colour = colour_set[22], method = "loess",
        span = 0.06, size = 0.475, se = FALSE
      ) +
      geom_point(
        aes(pos_mb, base, colour = type), size = 1.5, shape = 25
      ) +
      scale_colour_manual(
        legend_title, values = colours_comps_genes,
        limits = levels(as.factor(snp_data$type)),
        guide = guide_legend(
          override.aes = list(
            shape = c(rep(16, 4), 25, 25)
          )
        )
      )
  }
)

# turn plot list into ggmatrix
plots_matrix <- ggmatrix(
  plots,
  nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = str_c(
    "Normalized Jost's D Value"
  ),
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = "Marker-by-Marker Jost's D Coloured by Distinct Group",
  legend = c(1, 1)
)

# plot the matrix
png(
  file.path("results", "phys_pos_vs_josts_D_with_distinct_group.png"),
  family = "Times New Roman", width = 320, height = 240, units = "mm", res = 192
)
plots_matrix + theme(
  legend.position = "bottom", legend.box = "vertical",
  plot.title = element_text(size = 20), axis.title = element_text(size = 20),
  legend.text = element_text(size = 16)
)
dev.off()

################################################################################
# print out all markers on each chromosome
blah <- by(
  snp_data %>% select(id, chrom, pos_mb, chrs, chrw, csws, josts_d, type),
  snp_data$chrom,
  function(chrom) {
    write_csv(
      chrom,
      file.path(josts_d_by_chrom,
        str_c(
          chrom$chrom[1],
          "_all_markers_Ds_and_major_allele_freqs_by_pop_with_genes.csv"
        )
      )
    )
  }
)
