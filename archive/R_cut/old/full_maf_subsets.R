library(SNPRelate)

setwd("C:\\Users\\Max_H\\OneDrive - University of Guelph\\Pedagogy\\PGDA\\")

gdsSubset <- "Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds"
source("Analysis\\R\\functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\wheat_phys_subset_both.gds")

## remove the minor alleles
for (maf in c(seq(0, 0.2, 0.01), seq(.25, .45, .05))) {
  set.seed(1000)
  snpset1 <- snpgdsSelectSNP(wheat, autosome.only = F, maf = maf)
  set.seed(1000)
  snpset2 <- snpgdsSelectSNP(wheat, autosome.only = F, maf = maf + 0.05)
  snp.indices1 <- match(unlist(snpset1), snp.id)
  snp.indices2 <- match(unlist(snpset2), snp.id)
  snp.indices <- snp.indices1[!snp.indices1 %in% snp.indices2]
  
  snpgdsCreateGeno(paste("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_slice_", maf*100, "_to_", (maf+0.05)*100, ".gds", sep = ""),
                   genmat = genotypes[snp.indices,],
                   sample.id = sample.id,
                   snp.id = snp.id[snp.indices],
                   snp.chromosome = snp.chrom[snp.indices],
                   snp.position = snp.pos[snp.indices],
                   other.vars = list(samp.annot = samp.annot),
                   snpfirstdim = T)
  snpgdsCreateGeno(paste("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_minus_", maf*100 ,".gds", sep = ""),
                   genmat = genotypes[-snp.indices1,],
                   sample.id = sample.id,
                   snp.id = snp.id[-snp.indices1],
                   snp.chromosome = snp.chrom[-snp.indices1],
                   snp.position = snp.pos[-snp.indices1],
                   other.vars = list(samp.annot = samp.annot),
                   snpfirstdim = T)
  snpgdsCreateGeno(paste("Data\\Intermediate\\GDS\\wheat_phys_subset_both_maf_", maf*100 ,".gds", sep = ""),
                   genmat = genotypes[snp.indices1,],
                   sample.id = sample.id,
                   snp.id = snp.id[snp.indices1],
                   snp.chromosome = snp.chrom[snp.indices1],
                   snp.position = snp.pos[snp.indices1],
                   other.vars = list(samp.annot = samp.annot),
                   snpfirstdim = T)
}
snpgdsClose(wheat)
