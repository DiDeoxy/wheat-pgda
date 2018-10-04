library(tidyverse)
library(SNPRelate)

parse_gds <- function(gds_name) {
  gds <- snpgdsOpen(str_c("Data\\Intermediate\\GDS\\", gds_name, ".gds"))
  # make a tibble of the snp data
  snp <- tibble(
    id = as.character(read.gdsn(index.gdsn(gds, "snp.id"))),
    chrom = as.integer(
      read.gdsn(index.gdsn(gds, "snp.chromosome"))
    ),
    pos = as.integer(read.gdsn(index.gdsn(gds, "snp.position"))),
    pos_mb = as.integer(read.gdsn(index.gdsn(gds, "snp.position"))) / 1e6
  )

  # extract the genotypes
  genotypes <- read.gdsn(index.gdsn(gds, "genotype"))

  # extract the categorical information
  sample_annot <- read.gdsn(index.gdsn(gds, "samp_annot"))
  sample_id <- as.character(read.gdsn(index.gdsn(gds, "sample.id")))
  snpgdsClose(gds)

  return(
    list(
      snp = snp, genotypes = genotypes,
      sample = list(id = sample_id, annot = sample_annot)
    )
  )
}

snpgds_create_snp_subset <- function(wheat_data, subset, snp_index) {
  snpgdsCreateGeno(
    str_c("Data\\Intermediate\\GDS\\", subset, ".gds"),
    genmat = wheat_data$genotypes[snp_index, ],
    sample.id = wheat_data$sample$id,
    snp.id = wheat_data$snp$id[snp_index],
    snp.chromosome = wheat_data$snp$chrom[snp_index],
    snp.position = wheat_data$snp$pos[snp_index],
    other.vars = list(samp_annot = wheat_data$sample$annot),
    snpfirstdim = T
  )
}

snpgds_create_sample_subset <- function(wheat_data, subset, sample_index) {
  snpgdsCreateGeno(
    str_c("Data\\Intermediate\\GDS\\", subset, ".gds"),
    genmat = wheat_data$genotypes[, -sample_index],
    sample.id = wheat_data$sample$id[-sample_index],
    snp.id = wheat_data$snp$id,
    snp.chromosome = wheat_data$snp$chrom,
    snp.position = wheat_data$snp$pos,
    other.vars = list(samp_annot = wheat_data$sample$annot),
    snpfirstdim = T
  )
}
