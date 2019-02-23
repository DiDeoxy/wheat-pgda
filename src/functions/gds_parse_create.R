library(tidyverse)
library(SNPRelate)

base <- file.path("data", "R", "GDS")

parse_gds <- function(gds_file) {
  gds <- snpgdsOpen(gds_file)

  # make a tibble of the snp data
  snp <- tibble(
    id = as.character(read.gdsn(index.gdsn(gds, "snp.id"))),
    chrom = read.gdsn(index.gdsn(gds, "snp.chromosome")),
    pos = read.gdsn(index.gdsn(gds, "snp.position")),
    pos_mb = read.gdsn(index.gdsn(gds, "snp.position")) / 1e6
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
      snp = snp, genotypes = genotypes,
      sample = list(id = sample_id, annot = sample_annot)
    )
  )
}

snpgds_create_snp_subset <- function(wheat_data, gds_file, snp_index) {
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

snpgds_create_sample_subset <- function(wheat_data, gds_file, sample_index) {
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
