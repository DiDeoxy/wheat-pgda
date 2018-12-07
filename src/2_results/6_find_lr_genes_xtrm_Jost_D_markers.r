library(tidyverse)

# load custom functions
source("src/R_functions/funcs_gds_parse_create.R")
source("src/R_functions/funcs_locus_by_locus.r")

# load the data from the gds object
wheat_data <- parse_gds("mr_pruned_phys_sample_subset")

# find the max position of any marker on each genome for xlims
max_genome_lengths <- calc_max_genome_lengths(wheat_data)

# names of groups to be plotted
groups <- c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")

# find the phi values of each marker in each Gene and add to data set
wheat_data <- add_group_stat(wheat_data, groups)

# find the extremeicance threshold for each Gene
extremes <- calc_extremes(wheat_data, groups)

group_extreme_freqs <- calc_group_extreme_freqs(
  wheat_data, extremes, prune = TRUE
)

# load the gene positions
known_genes <- load_groups("known_genes.csv")[c(1, 6, 9, 11, 15), ]
# add min phi values of plotting to each gene so they appear at bottom of plots
known_genes <- cbind(known_genes, base = 0)
  
base <- "Results/loci/D/closest_markers"

for (i in 1:nrow(known_genes)) {
  row <- known_genes[i, ]
  ifelse(! dir.exists(file.path(base, row$chrom[1])),
    dir.create(file.path(base, row$chrom[1])), FALSE
  )
  chrom_group_extreme_freqs <- group_extreme_freqs %>%
    filter(chrom == row$chrom & group == row$group) %>%
    drop_na(freq) %>%
    bind_rows(row) %>%
    arrange(pos_mb)
  write_csv(chrom_group_extreme_freqs,
    file.path(base, row$chrom, str_c(row$gene, "_high_freq_extreme.csv"))
  )
}