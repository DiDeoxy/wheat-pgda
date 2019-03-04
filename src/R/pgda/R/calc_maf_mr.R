#' Calculate the MAF and MR
#'
#' Calculates the minor allele frequency (MAF) and missing rate (MR) of each
#' marker
#'
#' @importFrom dplyr bind_rows tibble
#' @importFrom magrittr %>%
#'
#' @param wheat_data the parsed gds object
#'
#' @return a list contianing the above data organised by genome

calc_maf_mr <- function (wheat_data) {
  chrom_geno_sums <- by(wheat_data$genotypes, wheat_data$snp$chrom,
    function (chrom_genos) {
      apply(chrom_genos, 1, function (snp) {
        A <- sum(snp == 0)
        B <- sum(snp == 2)
        missing <- sum(snp == 3)
        return(tibble(maf = min(c(A, B)) / sum(A, B), mr = missing / sum(A, B, missing)))
      }) %>% bind_rows()
    }
  )
  list(
    A = chrom_geno_sums[seq(1, 19, 3)] %>% bind_rows(),
    B = chrom_geno_sums[seq(2, 20, 3)] %>% bind_rows(),
    D = chrom_geno_sums[seq(3, 21, 3)] %>% bind_rows()
  )
}