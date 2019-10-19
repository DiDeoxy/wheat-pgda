library(tidyverse)
library(poppr)
library(SNPRelate)
library(ape)

## loading the gds of the data and pullling some attributes out
source("src\\R_functions\\funcs_gds_parse_create.R")

wheat_data <- parse_gds("phys_subset_sample")

wheat_gds <- snpgdsOpen(
  "Data\\Intermediate\\GDS\\phys_subset_sample_ld_pruned.gds"
)
sample_dist <- 1 - snpgdsIBS(wheat_gds, autosome.only = F)$ibs
snpgdsClose(wheat_gds)

comparisons <- c(
  "chrs_csws", "chrs_chrw", "csws_chrw", "Lr1", "Lr10", "Lr21", "Lr22a", "Lr34"
)

for (comparison in comparisons) {
  genind <- read_rds(
    str_c("Data\\Intermediate\\Adegenet\\", comparison, "_genind.rds")
  )

  genind_seploc <- seploc(genind)
  for (i in 1:length(genind_seploc)) {
    strata(genind_seploc[[i]]) <- strata(genind)
  }
  subset_amovas <- list()
  for (i in 1:length(genind_seploc)) {
    if (isPoly(genind_seploc[[i]])) {
      locus <- missingno(genind_seploc[[i]], type = "geno")
      indivs <- match(rownames(locus$tab), wheat_data$sample$id)
      subset_amovas[[i]] <- poppr.amova(
        locus,
        hier = str_c("~", comparison) %>% as.formula(), missing = "genotype",
        within = F, clonecorrect = F,
        dist = as.dist(sample_dist[indivs, indivs])
      )
    } else {
      subset_amovas[[i]] <- NA
    }
  }

  write_rds(subset_amovas,
    path = str_c("Data\\Intermediate\\Adegenet\\", comparison, "_amova.rds")
  )
}