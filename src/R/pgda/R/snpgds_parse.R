#' Parse a gds file
#'
#' Loads the data from a gds object into a form amenable to manipulation
#' using R code
#'
#' @importFrom dplyr tibble
#' @importFrom magrittr %>%
#' @importFrom SNPRelate snpgdsClose snpgdsOpen
#' @importFrom gdsfmt index.gdsn read.gdsn
#'
#' @param gds_file the path to the target gds file
#'
#' @return a list of lists containing the data of the gds object
#'
#' @export

snpgds_parse <- function(gds_file) {
  gds <- snpgdsOpen(gds_file)

  # make a tibble of the snp data
  snp <- tibble(
    id = read.gdsn(index.gdsn(gds, "snp.id")) %>% as.character(),
    chrom = read.gdsn(index.gdsn(gds, "snp.chromosome")),
    pos = read.gdsn(index.gdsn(gds, "snp.position")),
    pos_mb = read.gdsn(index.gdsn(gds, "snp.position")) / 1e6
  )

  max_lengths <- by(snp$pos, snp$chrom, max) %>%
    (
      function (max_chrom_lengths) {
        c(A = max(max_chrom_lengths[seq(1, 19, 3)]),
          B = max(max_chrom_lengths[seq(2, 20, 3)]),
          D = max(max_chrom_lengths[seq(3, 21, 3)]), 
          one = max(max_chrom_lengths[c(1, 2, 3)]),
          two = max(max_chrom_lengths[c(4, 5, 6)]),
          three = max(max_chrom_lengths[c(7, 8, 9)]),
          four = max(max_chrom_lengths[c(10, 11, 12)]),
          five = max(max_chrom_lengths[c(13, 14, 15)]),
          six = max(max_chrom_lengths[c(16, 17, 18)]),
          seven = max(max_chrom_lengths[c(19, 20, 21)])
        )
      }
    )

  # convert chroms values from integer to names
  chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>% 
    t() %>% as.vector()
  for (i in 1:21) {
    snp$chrom[which(snp$chrom == i)] <- chroms[i]
  }

  # extract the genotypes
  genotypes <- read.gdsn(index.gdsn(gds, "genotype"))

  # extract the categorical information
  sample_annot <- read.gdsn(index.gdsn(gds, "samp_annot"))
  sample_id <- as.character(read.gdsn(index.gdsn(gds, "sample.id")))
  snpgdsClose(gds)

  return(
    list(
      snp = snp, genotypes = genotypes, max_lengths = max_lengths,
      sample = list(id = sample_id, annot = sample_annot)
    )
  )
}