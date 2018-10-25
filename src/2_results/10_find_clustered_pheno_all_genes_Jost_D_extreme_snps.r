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

# add columns of the D values of each comparison
for (group in c("chrs_chrw", "chrs_csws", "csws_chrw")) {
  comp_Jost_D <- read_rds(
    str_c("Data\\Intermediate\\mmod\\", group, "_Jost_D.rds")
  )[[1]]
  wheat_data$snp <- wheat_data$snp %>% add_column(!!group := comp_Jost_D)
}

# find the top 2.5% quantile for each of the three comparisons
extremes <- wheat_data$snp %>%
  summarise(
    chrs_csws = quantile(chrs_csws, prob = 0.975, na.rm = T),
    chrs_chrw = quantile(chrs_chrw, prob = 0.975, na.rm = T),
    csws_chrw = quantile(csws_chrw, prob = 0.975, na.rm = T)
  )

# load gene data
pheno_genes <- read_csv(
  "Data\\Intermediate\\Aligned_genes\\selected_alignments\\pheno_genes.csv",
  col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
) %>%
  select(id, chrom, pos) %>%
  mutate(pos_mb = pos / 1e6, comparison = "pheno_gene")
resi_genes <- read_csv(
  "Data\\Intermediate\\Aligned_genes\\selected_alignments\\resi_genes.csv",
  col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
) %>%
  select(id, chrom, pos) %>%
  mutate(pos_mb = pos / 1e6, comparison = "resi_gene")
# join them together
genes <- pheno_genes %>% full_join(resi_genes)

# add min phi values of plotting to each gene so they appear at bottom of plots
genes <- add_column(
  genes, chrs_csws = 0, chrs_chrw = 0, csws_chrw = 0, dist = 20
)

base <- "Results\\loci\\D\\closest_markers"

for (i in 1:nrow(genes)) {
  row <- genes[i, ]
  nearby_markers <- wheat_data$snp %>%
    filter(chrom == row$chrom & abs(pos_mb - row$pos_mb) < row$dist) %>%
    transmute(Marker = id, Chrom = chrom, Distance = pos_mb - row$pos_mb,
      chrs_csws, chrs_chrw, csws_chrw
    )
  nearby_extreme_markers <- nearby_markers
  nearby_extreme_markers$chrs_csws[
    which(nearby_extreme_markers$chrs_csws < extremes$chrs_csws)
  ] <- NA
  nearby_extreme_markers$chrs_chrw[
    which(nearby_extreme_markers$chrs_chrw < extremes$chrs_chrw)
  ] <- NA
  nearby_extreme_markers$csws_chrw[
    which(nearby_extreme_markers$csws_chrw < extremes$csws_chrw)
  ] <- NA
  nearby_extreme_markers <- nearby_extreme_markers[
    ! is.na(nearby_extreme_markers$chrs_csws) |
    ! is.na(nearby_extreme_markers$chrs_chrw) |
    ! is.na(nearby_extreme_markers$csws_chrw),
  ]
  ifelse(! dir.exists(file.path(base, row$chrom)),
    dir.create(file.path(base, row$chrom)), FALSE)
  write_csv(nearby_markers,
    file.path(base, row$chrom, str_c(row$id, "_comparisons_full.csv"))
  )
  write_csv(nearby_extreme_markers,
    file.path(base, row$chrom, str_c(row$id, "_comparisons_extreme.csv"))
  )
}

# output all extreme markers sliced in various combinations
wheat_data$snp$chrs_csws[wheat_data$snp$chrs_csws < extremes$chrs_csws] <- NA
wheat_data$snp$chrs_chrw[wheat_data$snp$chrs_chrw < extremes$chrs_chrw] <- NA
wheat_data$snp$csws_chrw[wheat_data$snp$csws_chrw < extremes$csws_chrw] <- NA

wheat_data$snp <- wheat_data$snp[
  ! is.na(wheat_data$snp$chrs_csws) |
  ! is.na(wheat_data$snp$chrs_chrw) |
  ! is.na(wheat_data$snp$csws_chrw),
]

blah <- by(wheat_data$snp[, -3], wheat_data$snp$chrom, function (chrom) {
  ifelse(! dir.exists(file.path(base, chrom$chrom[1])),
    dir.create(file.path(base, chrom$chrom[1])), FALSE
  )
  ifelse(
    nrow(genes[genes$chrom == chrom$chrom[1], ]),
    chrom <- chrom %>%
      full_join(
        genes[genes$chrom == chrom$chrom[1], ]) %>%
      arrange(pos_mb) %>%
      select(-c(comparison, dist, pos)),
    FALSE
  )
  # one at a time
  write_csv(
    chrom[
      which(
        (! is.na(chrom$chrs_csws) | chrom$chrs_csws == 0) &
        (is.na(chrom$chrs_chrw) | chrom$chrs_chrw == 0) &
        (is.na(chrom$csws_chrw) | chrom$csws_chrw == 0)
      ),
    ],
    file.path(base, chrom$chrom[1], "chrs_vs_csws_extreme.csv")
  )
  write_csv(
    chrom[
      which(
        (is.na(chrom$chrs_csws) | chrom$chrs_csws == 0) &
        (! is.na(chrom$chrs_chrw) | chrom$chrs_chrw == 0) &
        (is.na(chrom$csws_chrw) | chrom$csws_chrw == 0)
      ),
    ],
    file.path(base, chrom$chrom[1], "chrs_vs_chrw_extreme.csv")
  )
  write_csv(
    chrom[
      which(
        (is.na(chrom$chrs_csws) | chrom$chrs_csws == 0) &
        (is.na(chrom$chrs_chrw) | chrom$chrs_chrw == 0) &
        (! is.na(chrom$csws_chrw) | chrom$csws_chrw == 0)
      ),
    ],
    file.path(base, chrom$chrom[1], "csws_vs_chrw_extreme.csv")
  )
  # pairs
  write_csv(
    chrom[
      which(
        (! is.na(chrom$chrs_chrw) | chrom$chrs_chrw == 0) &
        (! is.na(chrom$csws_chrw) | chrom$csws_chrw == 0)
      ),
    ],
    file.path(base, chrom$chrom[1],
      "chrs_vs_chrw_and_csws_vs_chrw_extreme.csv")
  )
  write_csv(
    chrom[
      which(
        (! is.na(chrom$chrs_csws) | chrom$chrs_csws == 0) &
        (! is.na(chrom$chrs_chrw) | chrom$chrs_chrw == 0) 
      ),
    ],
    file.path(base, chrom$chrom[1],
      "chrs_vs_csws_and_chrs_vs_chrw_extreme.csv")
  )
  write_csv(
    chrom[
      which(
        (! is.na(chrom$chrs_csws) | chrom$chrs_csws == 0) &
        (! is.na(chrom$csws_chrw) | chrom$csws_chrw == 0)
      ),
    ],
    file.path(base, chrom$chrom[1],
      "chrs_vs_csws_and_csws_vs_chrw_extreme.csv")
  )
})
