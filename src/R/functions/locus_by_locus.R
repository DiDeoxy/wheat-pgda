load_genes <- function(csv, base = 0) {
  read_csv(
    csv, col_names = c("id", "chrom", "pos", "sleng", "salign", "%id")
  ) %>%
    mutate(pos_mb = pos / 1e6) %>%
    select(id, chrom, pos_mb) %>%
    cbind(base = base)
}

calc_eh <- function (genotypes) {
  apply(genotypes, 1, function (snp) {
    2 * ((sum(snp == 0) / sum(snp == 0 | 2)) *
      (sum(snp == 2) / sum(snp == 0 | 2)))
  })
}