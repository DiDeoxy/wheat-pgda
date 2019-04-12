# load file paths, colours, and functions
source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "arrange", "mutate", "rowwise", "select")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_point", "geom_smooth", "guide_legend",
  "scale_colour_manual", "theme", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(pgda, "load_genes", "snpgds_parse")
import::from(plyr, "rbind.fill")
import::from(readr, "read_rds", "write_csv")
import::from(stringr, "str_c")
import::from(tibble, "add_column", "add_row", "as_tibble", "tibble")
import::from(tidyr, "gather")

# load the data from the gds object
wheat_data <- snpgds_parse(phys_gds)

# get the clusters
cluster <- read_rds(hdbscan)$cluster

index_chrs_chrw_csws <- which(
  (wheat_data$sample$annot$pheno == "HRS" & cluster == 5)
  | (wheat_data$sample$annot$pheno == "HRW" & cluster == 1)
  | (wheat_data$sample$annot$pheno == "SWS" & cluster == 2)
)
cpheno <- wheat_data$sample$annot$pheno[index_chrs_chrw_csws] %>% factor()
geno <- wheat_data$genotype[, index_chrs_chrw_csws]

major_allele_freq_pops <- apply(geno, 1, function(marker) {
  total_geno <- table(marker)
  major <- which.max(total_geno) %>% names()
  by(marker, cpheno, function(pop) {
    pop_geno <- table(pop)
    if (all(c("0", "2") %in% names(pop_geno))) {
      pop_geno[[major]] / sum(pop_geno[["0"]], pop_geno[["2"]])
    } else if (major %in% names(pop_geno)) {
      1
    } else {
      0
    }
  }) %>% rbind()
}) %>% t() %>% as_tibble()
colnames(major_allele_freq_pops) <- c("chrs", "chrw", "csws")

group <- "chrs_chrw_csws"
# add the genes positons to the regions table
wheat_data$snp <- wheat_data$snp %>%
  add_column(group := read_rds(josts_d)) %>%
  gather(group, D, group) %>%
  cbind(major_allele_freq_pops %>% round(4)) %>%
  rowwise() %>%
  mutate(class =
    ifelse(
      D > 0.32,
      c("CHRSD", "CHRWD", "CSWSD")[
        which.max(
          c(
            sum(abs(chrs - chrw), abs(chrs - csws)),
            sum(abs(chrw - chrs), abs(chrw - csws)), 
            sum(abs(csws - chrs), abs(csws - chrw))
          )
        )
      ],
      "None"
    )
  )
    

# overall freqs
summary(wheat_data$snp$D)
sum(wheat_data$snp$D > mean(wheat_data$snp$D)) / length(wheat_data$snp$D)

# genome freqs
for (genome in c("A", "B", "D")) {
  median(wheat_data$snp$D[which(grepl(genome, wheat_data$snp$chrom))]) %>%
    round(4) %>% str_c(genome, " median = ", .) %>% print()
}

# chromsome group freqs
for (chr_group in 1:7) {
  median(wheat_data$snp$D[which(grepl(chr_group, wheat_data$snp$chrom))]) %>%
    round(4) %>% str_c("Chr ", chr_group, " median = ", .) %>% print()
}

# group freq and median values
classes <- c("CHRSD", "CHRWD", "CSWSD", "None")

for (class in classes) {
  (sum(wheat_data$snp$class == class) / length(wheat_data$snp$class) * 100) %>%
    round(2) %>% str_c("Class ", class, " percent = ", .) %>% print()
  median(wheat_data$snp$D[which(wheat_data$snp$class == class)]) %>%
    round(2) %>% str_c("Class ", class, " median = ", .) %>% print()
}

# load the gene positions
pheno_genes <- load_genes(
  file.path(blast, "selected_pheno.csv"), base = 1
) %>% mutate(gene_type = "Phenotype Genes")
resi_genes <- load_genes(
  file.path(blast, "selected_resi.csv"), base = 1
) %>% mutate(gene_type = "Resistance Genes")

wheat_data$snp <- wheat_data$snp %>%
  rbind.fill(pheno_genes, resi_genes) %>%
  arrange(chrom, pos_mb) %>%
  mutate(type = pmin(gene_type, class, na.rm = TRUE)) %>%
  select(-c(gene_type, class))

# create a list of plots, one for each chromosome with the correct markers and
# genes on each coloured by comparison or gene type
legend_title <- "Marker Types & Genes"
plots <- by(
  wheat_data$snp, wheat_data$snp$chrom, function(chrom_groups) {
    chrom <- chrom_groups$chrom[1]
    chrom_groups %>%
      ggplot() +
      ylim(0, 1) +
      xlim(
        0,
        wheat_data$max_lengths[
          ifelse(grepl("A", chrom), 1, ifelse(grepl("B", chrom), 2, 3))
        ] / 1e6
      ) +
      geom_point(
        aes(pos_mb, D, colour = type), shape = 16, size = 1
      ) +
      geom_smooth(
        aes(pos_mb, D), colour = colour_set[22], method = "loess", span = 0.06,
        size = 0.475, se = FALSE
      ) +
      geom_point(
        aes(pos_mb, base, colour = type), size = 1.5, shape = 25
      ) +
      scale_colour_manual(
        legend_title, values = colours_comps_genes,
        limits = levels(as.factor(wheat_data$snp$type)),
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
  title = "Normalized Jost's D Values with Loess Curve By Chromosome ",
  legend = c(1, 1)
)

# plot the matrix
png(
  file.path("results", "clustered_phenos_markers_josts_ds_with_genes.png"),
  family = "Times New Roman", width = 210, height = 267, pointsize = 5,
  units = "mm", res = 300
)
plots_matrix + theme(legend.position = "bottom", legend.box = "vertical")
dev.off()

################################################################################
# print out all markers on each chromosome
blah <- by(wheat_data$snp %>% select(-c(base, group)), wheat_data$snp$chrom,
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