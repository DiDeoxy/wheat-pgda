library(tidyverse)
library(SNPRelate)

# load data into R object from the GDS object using the script at the source
# location
gds <- "Data\\Intermediate\\GDS\\full_phys.gds"
source("src\\functions\\data_loading.R")

# open the GDS object again and select SNPs with a missing data rate below 0.1
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys.gds")
snpset_ids_list <- snpgdsSelectSNP(
    wheat, autosome.only = F, missing.rate = 0.1)
snpgdsClose(wheat)
snp_indices <- match(unlist(snpset_ids_list), snp_id)

# create a new gds object with only the snps with less than 10% missing data
snpgdsCreateGeno("Data\\Intermediate\\GDS\\full_phys_snpgdsSelectSNP.gds",
                 genmat = genotypes[snp_indices,],
                 sample.id = sample_id, 
                 snp.id = snp_id[snp_indices],
                 snp.chromosome = snp_chrom[snp_indices],
                 snp.position = snp_pos[snp_indices],
                 other.vars = list(samp_annot = samp_annot),
                 snpfirstdim = T)

# same aas above but for gen map gds
gds <- "Data\\Intermediate\\GDS\\full_gen.gds"
source("src\\functions\\data_loading.R")
snpgdsCreateGeno("Data\\Intermediate\\GDS\\full_gen_snpgdsSelectSNP.gds",
                 genmat = genotypes[snp_indices,],
                 sample.id = sample_id, 
                 snp.id = snp_id[snp_indices],
                 snp.chromosome = snp_chrom[snp_indices],
                 snp.position = snp_pos[snp_indices],
                 other.vars = list(samp_annot = samp_annot),
                 snpfirstdim = T)

# eliminate those individuals that show identity by state 
# (IBS, fractional identity) greater than 0.99
gds <- "Data\\Intermediate\\GDS\\full_phys_snpgdsSelectSNP.gds"
source("src\\functions\\data_loading.R")

wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_snpgdsSelectSNP.gds")
IBS <- snpgdsIBS(wheat, autosome.only = F)
snpgdsClose(wheat)

pairs <- which(IBS$ibs >= 0.99, arr.ind = T)

indices <- vector()
for (i in 1:dim(pairs)[1]) {
  if (pairs[i,1] == pairs[i,2]){
    indices <<- c(indices, i)
  }
}
pairs <- pairs[-indices,]
pairs

# i take the pairs and use a perl function to find the connected graphs
# i.e. those sets of pairs that form a connected triangle, square or more

# eliminate individuals so that only one from each pair, triangle, or square
# remains
NILs <- c("PT434", "BW811", "AC Minto 1", "Avocet 1", "BW275 1", "BW395",
          "PT616", "BW427 1", "BW492", "BW948", "Carberry 1", "CDC Stanley 1",
          "PT754", "SWS349", "Somerset 1", "Stettler 1", "AC Reed 1", "SWS87",
          "SWS241", "SWS345", "SWS363", "SWS390", "SWS408", "SWS410")

# find the indices of the NILs
sample_indices <- match(NILs, sample_id)

# eliminate the nils from the sample annotation info
retained_samp_annot <-  list()
for(name in names(samp_annot)) {
  retained_samp_annot[[name]] <- factor(samp_annot[[name]][-sample_indices])
}
#create gds ubjects without the NILs, physical map
snpgdsCreateGeno("Data\\Intermediate\\GDS\\full_phys_subset_sample.gds",
                 genmat = genotypes[,-sample_indices],
                 sample.id = sample_id[-sample_indices], 
                 snp.id = snp_id,
                 snp.chromosome = snp_chrom,
                 snp.position = snp_pos,
                 other.vars = list(samp_annot = retained_samp_annot),
                 snpfirstdim = T)

# gen map
gds <- "Data\\Intermediate\\GDS\\full_gen_snpgdsSelectSNP.gds"
source("src\\functions\\data_loading.R")
snpgdsCreateGeno("Data\\Intermediate\\GDS\\full_gen_subset_sample.gds",
                 genmat = genotypes[,-sample_indices],
                 sample.id = sample_id[-sample_indices], 
                 snp.id = snp_id,
                 snp.chromosome = snp_chrom,
                 snp.position = snp_pos,
                 other.vars = list(samp_annot = retained_samp_annot),
                 snpfirstdim = T)

# identify and remove from the overall dataset LD pruned markers
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\functions\\data_loading.R")

# performs LD pruning to produce a set of SNPs that maximally represent the diversity
# of the genome with as little redundant info as possible used for estimation of relationships
# genome wide (i.e. clustering)
wheat <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_subset_sample.gds")
set.seed(1000)
snpset_ids_list <- snpgdsLDpruning(
    wheat, autosome.only = F, missing.rate = 0.1, maf = 0.05,
    ld.threshold = 0.6, slide.max.bp = 1e7)
snpgdsClose(wheat)
snp_indices <- match(unlist(snpset_ids_list), snp_id)

# create a subset gds object containg only the SNPs remaing after ld pruning
snpgdsCreateGeno("Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds",
                 genmat = genotypes[snp_indices,],
                 sample.id = sample_id, 
                 snp.id = snp_id[snp_indices],
                 snp.chromosome = snp_chrom[snp_indices],
                 snp.position = snp_pos[snp_indices],
                 other.vars = list(samp_annot = samp_annot),
                 snpfirstdim = T)