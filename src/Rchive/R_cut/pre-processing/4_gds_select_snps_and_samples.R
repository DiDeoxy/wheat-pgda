library(tidyverse)
library(SNPRelate)

source("src/R_functions/funcs_gds_parse_create.R")

# load data into R object from the GDS object using the script at the source
# location
wheat_data <- parse_gds("full_phys")

# open the GDS object again and select SNPs with a missing data rate below 0.1
wheat_gds <- snpgdsOpen("Data/Intermediate/GDS/full_phys.gds")
kept_id <- unlist(snpgdsSelectSNP(
  wheat_gds,
  autosome.only = F, missing.rate = 0.1
))
snpgdsClose(wheat_gds)
kept_index <- match(kept_id, wheat_data$snp$id)

# create a new gds object with only the snps with less than 10% missing data
snpgds_create_snp_subset(wheat_data, "phys_select_snp", kept_index)

# same aas above but for gen map gds
wheat_data <- parse_gds("full_gen")
kept_index <- sort(match(kept_id, wheat_data$snp$id))
snpgds_create_snp_subset(wheat_data, "gen_select_snp", kept_index)

# eliminate those individuals that show identity by state
# (IBS, fractional identity) greater than 0.99
wheat_gds <- snpgdsOpen(
  "Data/Intermediate/GDS/phys_select_snp.gds"
)
IBS <- snpgdsIBS(wheat_gds, autosome.only = F)
snpgdsClose(wheat_gds)

pairs <- which(IBS$ibs >= 0.99, arr.ind = T)

indices <- vector()
for (i in 1:dim(pairs)[1]) {
  if (pairs[i, 1] == pairs[i, 2]) {
    indices <- c(indices, i)
  }
}
pairs <- pairs[-indices, ]
pairs

# i take the pairs and use a perl function to find the connected graphs
# i.e. those sets of pairs that form a connected triangle, square or more

# eliminate individuals so that only one from each pair, triangle, or square
# remains
NILs <- c(
  "PT434", "BW811", "AC Minto 1", "Avocet 1", "BW275 1", "BW395",
  "PT616", "BW427 1", "BW492", "BW948", "Carberry 1", "CDC Stanley 1",
  "PT754", "SWS349", "Somerset 1", "Stettler 1", "AC Reed 1", "SWS87",
  "SWS241", "SWS345", "SWS363", "SWS390", "SWS408", "SWS410"
)


# find the indices of the NILs
sample_index <- match(NILs, wheat_data$sample$id)

# phys map
wheat_data <- parse_gds("phys_select_snp")
# eliminate the nils from the sample annotation info
for (name in names(wheat_data$sample$annot)) {
  wheat_data$sample$annot[[name]] <- factor(
    wheat_data$sample$annot[[name]][-sample_index]
  )
}
# create gds object without the NILs
snpgds_create_sample_subset(
  wheat_data, "phys_subset_sample", sample_index
)

# gen map
wheat_data <- parse_gds("gen_select_snp")
for (name in names(wheat_data$sample$annot)) {
  wheat_data$sample$annot[[name]] <- factor(
    wheat_data$sample$annot[[name]][-sample_index]
  )
}
snpgds_create_sample_subset(
  wheat_data, "gen_subset_sample", sample_index
)
