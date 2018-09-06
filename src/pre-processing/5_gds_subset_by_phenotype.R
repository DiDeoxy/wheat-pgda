library(SNPRelate)

gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\functions\\data_loading.R")

# creating gds subsets for major phenotypic types in the data set
# used in analyses of diversity within these groups
HRS <- which(retained_samp_annot$designation == "HRS")

load("Data\\Intermediate\\dbscan\\dbscan.RData")
dbscan$cluster[which(dbscan$cluster == 0)] <- dbscan$cluster[which(dbscan$cluster == 0)] + 6
dbscan6 <- factor(dbscan$cluster)

indexNotCanadian <- c(which(dbscan6 != 1))
HRS <- HRS[-which(HRS %in% indexNotCanadian) ]

HRSretained_samp_annot <- list()
for(name in names(samp_annot)) {
  HRSretained_samp_annot[[name]] <- factor(retained_samp_annot[[name]][HRS])
}

snpgdsCreateGeno("Data\\Intermediate\\GDS\\HRS_phys_subset_sample.gds",
                 genmat = genotypes[,HRS],
                 sample_id = sample_id[HRS],
                 snp_id = snp_id,
                 snp_chromosome = snp_chrom,
                 snp_position = snp_pos,
                 other.vars = list(samp_annot = HRSretained_samp_annot),
                 snpfirstdim = T)

HRW <- which(retained_samp_annot$designation == "HRW")
HRWretained_samp_annot <- list()
for(name in names(samp_annot)) {
  HRWretained_samp_annot[[name]] <- factor(retained_samp_annot[[name]][HRW])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\HRW_phys_subset_sample.gds",
                 genmat = genotypes[,HRW],
                 sample_id = sample_id[HRW],
                 snp_id = snp_id,
                 snp_chromosome = snp_chrom,
                 snp_position = snp_pos,
                 other.vars = list(samp_annot = HRWretained_samp_annot),
                 snpfirstdim = T)

SWS <- which(retained_samp_annot$designation == "SWS")
SWSretained_samp_annot <- list()
for(name in names(samp_annot)) {
  SWSretained_samp_annot[[name]] <- factor(retained_samp_annot[[name]][SWS])
}
snpgdsCreateGeno("Data\\Intermediate\\GDS\\SWS_phys_subset_sample.gds",
                 genmat = genotypes[,SWS],
                 sample_id = sample_id[SWS],
                 snp_id = snp_id,
                 snp_chromosome = snp_chrom,
                 snp_position = snp_pos,
                 other.vars = list(samp_annot = SWSretained_samp_annot),
                 snpfirstdim = T)