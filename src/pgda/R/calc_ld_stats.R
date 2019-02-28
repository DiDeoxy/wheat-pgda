#' Calculate composite ld between markers
#'
#' Calculates composite ld for the markers on each chromosome in a gds object,
#' returning a dist object and the LD between neighbouring markers
#'
#' @importFrom magrittr %>%
#' @importFrom pracma Diag
#' @import SNPRelate
#'
#' @param gds a gds object
#' @param snp_data the parsed snp data of the gds object
#'
#' @return a list contianing the dist object and neighbouring lds for each
#' chromosome ordered by genome

calc_ld_stats <- function (gds, snp_data) {
  wheat_internal <- snpgdsOpen(gds)
  # Calcualte ld between all snps on each chromosome
  ld_stats <- by(snp_data, snp_data$chrom, function (chrom) {
    ld_mat <- snpgdsLDMat(wheat_internal, method = "composite",
      slide = -1, snp.id = chrom$id)$LD %>%
      abs() 
    return(
      list(
        pw = ld_mat %>% as.dist() %>% as.vector(),
        nbs = ld_mat %>% Diag(., 1)
      )
    )
  })
  snpgdsClose(wheat_internal)

  genome_ld <- list(A = list(), B = list(), D = list())
  for (i in 1:length(ld_stats)) {
    if (i %in% seq(1, 19, 3)) {
      genome_ld$A$pw <- c(genome_ld$A$pw, ld_stats[[i]]$pw)
      genome_ld$A$nbs <- c(genome_ld$A$nbs, ld_stats[[i]]$nbs)
    } else if (i %in% seq(2, 20, 3)) {
      genome_ld$B$pw <- c(genome_ld$B$pw, ld_stats[[i]]$pw)
      genome_ld$B$nbs <- c(genome_ld$B$nbs, ld_stats[[i]]$nbs)
    } else {
      genome_ld$D$pw <- c(genome_ld$D$pw, ld_stats[[i]]$pw)
      genome_ld$D$nbs <- c(genome_ld$D$nbs, ld_stats[[i]]$nbs)
    }
  }

  genome_ld
}