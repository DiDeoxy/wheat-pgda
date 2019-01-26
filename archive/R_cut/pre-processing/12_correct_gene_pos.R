library(tidyverse)

corrected_gene_pos <- function(file_in, file_out) {
  genes <- read_csv(file_in,
    col_names = c(
      "id", "alignment", "chrom", "pos", "sleng",
      "aleng", "%id"
    )
  ) %>%
    select(id, chrom, pos) %>%
    arrange(chrom, pos)
  chroms <- as.vector(t(outer(as.character(1:7), c("A", "B", "D"), paste,
    sep = ""
  )))
  chroms_part2 <- str_c(chroms, "_part2")
  part2_start <- tibble(
    chrom = chroms_part2,
    start = c(
      471304005, 438720154, 452179604, 462376173,
      453218924, 462216879, 454103970, 448155269,
      476235359, 452555092, 451014251, 451004620,
      453230519, 451372872, 451901030, 452440856,
      452077197, 450509124, 450046986, 453822637,
      453812268
    )
  )
  genes <- genes %>%
    left_join(part2_start) %<>%
    mutate(pos = rowSums(cbind(.$pos, .$start), na.rm = TRUE) / 1e6) %<>%
    mutate(chrom = substr(.$chrom, 1, 2)) %<>%
    mutate(chrom = as.integer(factor(chrom, levels = chroms))) %>%
    select(id, chrom, pos)

  return(genes)
}

write_rds(
  corrected_gene_pos(
    "Data\\Intermediate\\Aligned_genes\\top_main_genes_selected.csv"
  ),
  "Data\\Intermediate\\Aligned_genes\\top_main_genes_selected_corrected.rds"
)
write_rds(
  corrected_gene_pos(
    "Data\\Intermediate\\Aligned_genes\\top_resistance_genes_selected.csv"
  ),
  str_c(
    "Data\\Intermediate\\Aligned_genes\\top_resistance_genes_selected",
    "_corrected.rds"
  )
)
write_rds(
  corrected_gene_pos(
    "Data\\Intermediate\\Aligned_genes\\known_genes_groups.csv"
  ),
  "Data\\Intermediate\\Aligned_genes\\known_genes_groups_corrected.rds"
)