#' Create a gds file that contains a subset of the snps of another gds file
#'
#' Takes the data parsed from a gds file using snpgds_parse and creates a new
#' gds file without the snps given by snpindex
#'
#' @import SNPRelate
#'
#' @param wheat_data the parsed gds data
#' @param gds_file the path and name of the new gds file
#' @param snp_index the indices of the snps to be removed
#'
#' @return None
#'
#' @export

snpgds_snp_subset <- function(wheat_data, gds_file, snp_index) {
  snpgdsCreateGeno(
    gds_file,
    genmat = wheat_data$genotypes[snp_index, ],
    sample.id = wheat_data$sample$id,
    snp.id = wheat_data$snp$id[snp_index],
    snp.chromosome = as.integer(as.factor(wheat_data$snp$chrom[snp_index])),
    snp.position = wheat_data$snp$pos[snp_index],
    other.vars = list(samp_annot = wheat_data$sample$annot),
    snpfirstdim = T
  )
}
