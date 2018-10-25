library(tidyverse)

# load custom functions
source("src\\R_functions\\funcs_gds_parse_create.R")

# load the data from the gds object
wheat_data <- parse_gds("phys_subset_sample")

chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>%
  t() %>% as.vector()
for (i in seq_along(chroms)) {
  wheat_data$snp$chrom[wheat_data$snp$chrom == i] <- chroms[i]
}

# find the max position of any marker on each genome for xlims
chrom_lengths <- by(wheat_data$snp$pos_mb, wheat_data$snp$chrom, max)
max_genome_lengths <- data.frame(
  A = max(chrom_lengths[seq(1, 19, 3)]),
  B = max(chrom_lengths[seq(2, 20, 3)]),
  D = max(chrom_lengths[seq(3, 21, 3)])
)

# find the phi values of each marker in each Gene and add to data set
for (group in c("Lr34", "Lr22a", "Lr21", "Lr10", "Lr1")) {
  gene_Jost_D <- read_rds(str_c("Data\\Intermediate\\mmod\\", group,
    "_Jost_D.rds"))[[1]]
  wheat_data$snp <- wheat_data$snp %>% add_column(!!group := gene_Jost_D)
}

# find the extremeicance threshold for each Gene
extremes <- wheat_data$snp %>%
  summarise(
    Lr34 = quantile(Lr34, prob = 0.975, na.rm = TRUE),
    Lr22a = quantile(Lr22a, prob = 0.975, na.rm = TRUE),
    Lr21 = quantile(Lr21, prob = 0.975, na.rm = TRUE),
    Lr10 = quantile(Lr10, prob = 0.975, na.rm = TRUE),
    Lr1 = quantile(Lr1, prob = 0.975, na.rm = TRUE)
  )

# load the gene positions
known_genes <- read_csv(
  "Data\\Intermediate\\Aligned_genes\\selected_alignments\\known_genes.csv",
  col_names = c("Gene", "chrom", "pos", "sleng", "salign", "%id")
) %>%
  select(Gene, chrom, pos) %>%
  mutate(pos_mb = pos / 1e6)

# add min phi values of plotting to each gene so they appear at bottom of plots
known_genes <- cbind(known_genes, min_extreme = min(extremes))
known_genes <- cbind(known_genes, gene_extreme = 0)
for (i in 1:nrow(known_genes)) {
  known_genes$gene_extreme[i] <- extremes[, known_genes$Gene[i]]
}

known_genes_responsible <- cbind(
  known_genes[c(1, 6, 9, 11, 15), ], dist = c(6.11, 6.11, 11.25, 1.64, 11.25)
) %>% as.tibble()

base <- "Results\\loci\\D\\closest_markers"

for (i in 1:nrow(known_genes_responsible)) {
  row <- known_genes_responsible[i, ]
  ifelse(! dir.exists(file.path(base, row$chrom[1])),
    dir.create(file.path(base, row$chrom[1])), FALSE
  )
  nearby_markers <- wheat_data$snp %>%
    filter(chrom == row$chrom & abs(pos_mb - row$pos_mb) < row$dist) %>%
    transmute(Marker = id, Chrom = chrom, Distance = pos_mb - row$pos_mb,
      jost_D := !!as.symbol(row$Gene))
  nearby_extreme_markers <- wheat_data$snp %>%
    filter(chrom == row$chrom & abs(pos_mb - row$pos_mb) < row$dist &
      !!as.symbol(row$Gene) > row$gene_extreme) %>%
    transmute(Marker = id, Chrom = chrom, Distance = pos_mb - row$pos_mb,
      jost_D := !!as.symbol(row$Gene))
  write_csv(nearby_markers,
    file.path(base, row$chrom, str_c(row$Gene, "_lr_genes_nearby.csv"))
  )
  write_csv(nearby_extreme_markers,
    file.path(base, row$chrom, str_c(row$Gene, "_lr_genes_extreme.csv"))
  )
}