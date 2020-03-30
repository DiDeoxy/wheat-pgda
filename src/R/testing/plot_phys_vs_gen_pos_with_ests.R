source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "arrange", "mutate", "select")
import::from(
  pgda, "calc_eh", "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom"
)
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_point", "geom_vline", "labs",
  "scale_colour_manual", "theme", "unit", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(plyr, "rbind.fill")
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
  gen_pos_cm = gen_data$snp$pos[match(phys_data$snp$id, gen_data$snp$id)] / 100
)

A4L13_ests <- load_genes(
  file.path(blast, "top2_4AL13.csv"), base = 1
) %>% mutate(est_type = "4AL13 ESTs", est_pos_mb = pos / 1e6) %>%
  select(-pos)
A4L5_ests <- load_genes(
  file.path(blast, "top2_4AL5.csv"), base = 1
) %>% mutate(est_type = "4AL5 ESTs", est_pos_mb = pos / 1e6) %>%
  select(-pos)
A4L4_ests <- load_genes(
  file.path(blast, "top2_4AL4.csv"), base = 1
) %>% mutate(est_type = "4AL4 ESTs", est_pos_mb = pos / 1e6) %>%
  select(-pos)

snp_data <- snp_data %>%
  rbind.fill(A4L13_ests, A4L5_ests, A4L4_ests) %>%
  arrange(chrom)

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
    chrom_data$est_type <- factor(
      chrom_data$est_type, levels = unique(snp_data$est_type)
    )
    plot <- chrom_data %>%
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
      geom_point(aes(phys_pos_mb, gen_pos_cm), colour = "black", size = 0.5)
    if (any(! is.na(chrom_data$est_pos_mb))) {
      plot <- plot +
        geom_vline(aes(xintercept = est_pos_mb, colour = est_type), size = 0.5)
    }
    plot +
      scale_colour_manual(name = "EST Type", values = colour_set[c(1, 6, 4)])
  }
)

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb", ylab = "Position in cM",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_c(
    "Marker Pseudo-Chromosomal vs Genetic Position by\n",
    "Chromosome Coloured by Major Allele Frequency"
  ),
  legend = c(5, 3)
)

# plot the matrix
png(
  file.path("results", "phys_vs_gen_pos_with_ests.png"),
  family = "Times New Roman", width = 200, height = 240, pointsize = 5,
  units = "mm", res = 192)
plots_matrix + theme(legend.position = "bottom")
dev.off()