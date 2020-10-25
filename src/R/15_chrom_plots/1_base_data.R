source("wheat-pgda/src/R/file_paths.R")
source("wheat-pgda/src/R/colour_sets.R")
import::from(dplyr, "bind_cols", "mutate", "rowwise", "select", "ungroup")
import::from(magrittr, "%>%", "%<>%")
import::from(pgda, "load_genes", "max_lengths", "snpgds_parse", "span_by_chrom")
import::from(plyr, "rbind.fill")
import::from(readr, "read_csv", "read_rds", "write_rds")
import::from(scrime, "rowTables")
import::from(stringr, "str_c")
import::from(tibble, "as_tibble", "enframe", "tibble")
import::from(tidyr, "gather")

################################################################################
# delineate chromosome names
chroms <- as.vector(
  t(outer(as.character(1:7), c("A", "B", "D"), paste, sep = ""))
)

# load the data from the gds object
phys_data <- snpgds_parse(phys_gds)
gen_data <- snpgds_parse(gen_gds)

################################################################################
# calc the lengths of the different genomes and homoeologous sets
max_phys_lengths <- span_by_chrom(
  phys_data$snp$chrom, phys_data$snp$pos
) %>% max_lengths() / 1e6
max_gen_lengths <- span_by_chrom(
  gen_data$snp$chrom, gen_data$snp$pos
) %>% max_lengths() / 100

################################################################################
# count genotypes
genos <- replace(phys_data$genotypes, phys_data$geno == 3, NA)

# calc allele counts
allele_counts <- rowTables(genos, c(0, 2))

# calc mjafs
mjafs <- do.call(pmax, (allele_counts / rowSums(allele_counts)) %>% as.data.frame())

# calc eh
eh <- (mjafs * (1 - mjafs)) * 2

# calc the mjafs by clustrer
cluster <- read_rds(hdbscan_rds)$cluster
cluster_genos <- list(
  CHRS = genos[,
    which(phys_data$sample$annot$mtg == "HRS" & cluster == 5)
  ],
  CHRW = genos[,
    which(phys_data$sample$annot$mtg == "HRW" & cluster == 1)
  ],
  CSWS = genos[,
    which(phys_data$sample$annot$mtg == "SWS" & cluster == 2)
  ]
)

mja <- apply(allele_counts, 1, which.max)
cluster_mjafs <- lapply(cluster_genos, function (sub_genos) {
  sub_allele_counts <- rowTables(sub_genos, c(0, 2))
  max_genos <- sub_allele_counts[cbind(seq_along(mja), mja)]
  max_genos / rowSums(sub_allele_counts)
}) %>% bind_cols()

################################################################################
# calc order diff inetervals
gen_to_phys_order <- match(phys_data$snp$id, gen_data$snp$id)
order_diffs <- (gen_to_phys_order - 1:length(gen_to_phys_order)) %>% abs()

order_diff_quantiles <- c("0%" = -1,
  quantile(
    order_diffs, c(1 / 2, 3 / 4, 7 / 8, 15 / 16, 31 / 32, 63 / 64, 127 / 128)
  ),
  "100%" = max(order_diffs)
)

intervals <- lapply(
  seq_along(order_diff_quantiles),
  function (i) {
    if (i < length(order_diff_quantiles)) {
      str_c(order_diff_quantiles[i] + 1, "-", order_diff_quantiles[i + 1])
    }
  }
) %>% unlist()

order_diff_intervals <- cut(order_diffs, order_diff_quantiles, intervals)

################################################################################
# collecting the data so far into one table
snp_data <- tibble(
  chrom = phys_data$snp$chrom,
  id = phys_data$snp$id,
  pos_mb = phys_data$snp$pos / 1e6,
  pos_cm = gen_data$snp$pos[gen_to_phys_order] / 100,
  eh = eh,
  mjaf = mjafs,
  odi = factor(order_diff_intervals),
  josts_d = read_rds(josts_d_rds)
) %>%
  bind_cols(cluster_mjafs) %>%
  rowwise() %>%
  mutate(
    josts_d_class = c("CHRSD", "CHRWD", "CSWSD")[
      which.max(
        c(
          sum(abs(CHRS - CHRW), abs(CHRS - CSWS)),
          sum(abs(CHRW - CHRS), abs(CHRW - CSWS)),
          sum(abs(CSWS - CHRS), abs(CSWS - CHRW))
        )
      )
    ]
  ) %>%
  ungroup() %>%
  mutate(josts_d_class = factor(josts_d_class))

# reorg mjafs and eh data for easier plotting
mjafs_by_cluster <- select(snp_data, chrom, pos_mb, mjaf, CHRS, CHRW, CSWS) %>%
  gather(cluster, mjaf, CHRS, CHRW, CSWS) %>%
  as_tibble() %>%
  split(.$chrom)

ehs_by_cluster <- select(snp_data, chrom, pos_mb, eh, CHRS, CHRW, CSWS) %>%
  gather(cluster, mjaf, CHRS, CHRW, CSWS) %>%
  as_tibble() %>%
  split(.$chrom)

snp_data %<>% split(.$chrom)

################################################################################
# create a tibble of the genes and centromeres

landmarks <- rbind.fill(
  load_genes(
    file.path(blast, "selected_pheno.csv")
  ) %>% mutate(
    type = "Gene", base = 0.75, pos_mb = pos / 1e6
  ) %>%
    dplyr::select(-pos),
  load_genes(
    file.path(blast, "selected_resi.csv")
  ) %>% mutate(
    type = "Gene", base = 0.75, pos_mb = pos / 1e6
  ) %>%
    dplyr::select(-pos),
  cbind(
    id = "Centromere", type = "Centromere", base = 0.75,
    read_csv(file.path(intermediate, "centromeres.csv"))
  )
) %>% as_tibble()

# or use go genes as landmarks
# landmarks <- read_csv(
#   file.path("/workspace", "data", "intermediate", "auxin_genes.csv")
# )

# landmarks <- read_csv(
#   file.path("/workspace", "data", "intermediate", "wounding_genes.csv")
# )