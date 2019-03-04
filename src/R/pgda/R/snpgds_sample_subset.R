#' Create a gds file that cotnains a subset of the sample from another gds file
#'
#' Takes the data parsed from a gds file using snpgds_parse and creates a new
#' gds file without the samples given by sample_index
#'
#' @import SNPRelate
#'
#' @param wheat_data the parsed gds data
#' @param gds_file the path and name of the new gds file
#' @param sample_index the indices of the samples to be removed
#'
#' @return None
#'
#' @export

snpgds_sample_subset <- function(wheat_data, gds_file, sample_index) {
  # eliminate the given samples from the annotations
  for (name in names(wheat_data$sample$annot)) {
    wheat_data$sample$annot[[name]] <- factor(
      wheat_data$sample$annot[[name]][-sample_index]
    )
  }
  
  snpgdsCreateGeno(
    gds_file,
    genmat = wheat_data$genotypes[, -sample_index],
    sample.id = wheat_data$sample$id[-sample_index],
    snp.id = wheat_data$snp$id,
    snp.chromosome = as.integer(as.factor(wheat_data$snp$chrom)),
    snp.position = wheat_data$snp$pos,
    other.vars = list(samp_annot = wheat_data$sample$annot),
    snpfirstdim = T
  )
}