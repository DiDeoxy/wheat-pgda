source(file.path("src", "R", "file_paths.R"))
source(file.path("src", "R", "colour_sets.R"))
import::from(dplyr, "select")
import::from(pgda, "calc_eh", "snpgds_parse")
import::from(GGally, "ggmatrix")
import::from(
  ggplot2, "aes", "ggplot", "geom_smooth", "labs", "scale_colour_manual",
  "theme", "unit", "xlim", "ylim"
)
import::from(magrittr, "%>%")
import::from(readr, "read_rds")
import::from(stringr, "str_c")
import::from(tibble, "tibble", "add_column")
import::from(tidyr, "gather")

# load the data from the gds object
wheat_data_phys <- snpgds_parse(phys_gds)
wheat_data_gen <- snpgds_parse(gen_gds)

snp_phys_order <- match(wheat_data_phys$snp$id, wheat_data_gen$snp$id)

clusters <- factor(read_rds(hdbscan)$cluster)
levels(clusters) <- c(
  "Noise", "Cluster 1 (HRW)", "Cluster 2 (SWS)", "Cluster 3 (CWES)",
  "Cluster 4 (CPSR/W)", "Cluster 5 (HRS)"
)

# make a tibble with the relevant data
phys_gen_snp_pos <- tibble(
  chrom = wheat_data_phys$snp$chrom, phys = wheat_data_phys$snp$pos_mb
)

for (cluster in levels(clusters)) {
  phys_gen_snp_pos <- phys_gen_snp_pos %>%
    add_column(
      !!cluster := calc_eh(
        wheat_data_phys$genotypes[, which(clusters == cluster)]
      )
    )
}
phys_gen_snp_pos <- phys_gen_snp_pos %>%
    add_column(All = calc_eh(wheat_data_phys$genotypes)) %>%
    select(-Noise) %>%
    gather(
      group, eh, `Cluster 1 (HRW)`, `Cluster 2 (SWS)`, `Cluster 3 (CWES)`,
      `Cluster 4 (CPSR/W)`, `Cluster 5 (HRS)`, All
    )

# allows application of same colour to each set of chromosomes
chroms_order <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>% 
  t() %>% as.vector()
colour_order <- c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3),
rep(6, 3), rep(7, 3))

# create a function for making a gradient of colours
colour_gradient <- colorRampPalette(c("Red", "Green", "Blue"))

plots <- by(phys_gen_snp_pos, phys_gen_snp_pos$chrom,
  function (data_chrom) {
    chrom <- data_chrom$chrom[1]
    data_chrom %>%
      ggplot() +
      xlim(
        0,
        wheat_data_phys$max_lengths[[
          ifelse(grepl("A", chrom), "A", ifelse(grepl("B", chrom), "B", "D"))
        ]] / 1e6
      ) +
      ylim(0, 0.5) +
      geom_smooth(aes(phys, eh, colour = group), size = 0.5, se = F) +
      scale_colour_manual(values = colours_hdbscan_legend) +
      labs(colour = "Cluster")
  }
)

plots_matrix <- ggmatrix(
  plots, nrow = 7, ncol = 3, xlab = "Position in Mb",
  ylab = "Expected Heterozygosity",
  xAxisLabels = c("A", "B", "D"), yAxisLabels = 1:7,
  title = str_c(
    "Comparions of physical position vs Expected\nHeterozygosity by cluster and All"
  ),
  legend = c(1, 1)
)

# plot the matrix
png(
  file.path("results", "physical_position_vs_eh_by_clusters_and_all.png"),
  family = "Times New Roman", width = 120, height = 267, pointsize = 5,
  units = "mm", res = 300)
plots_matrix + theme(legend.position = "bottom")
dev.off()