library(SNPRelate)

# load data
gds <- "Data\\Intermediate\\GDS\\full_phys.gds"
source("src\\R_functions\\data_loading.R")

source("src\\R_functions\\ld_floor_prune.R")

# create a list of the ld between markers on each chromosome add NA column
# and rows every million bps so that gaps appear on image
full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys.gds")
kept_ids <- by(data, snp_chrom, function(chrom) {
    ld_mat <- abs(snpgdsLDMat(full,
        method = "composite", snp.id = chrom$id,
        slide = -1
    )$LD)

    ld_floor <- 0.1
    max_window <- 0.1
    set.seed(1000)
    chrom$id[ld_floor_prune(chrom, ld_mat, ld_floor, max_window)]
})
snpgdsClose(wheat)
head(kept_ids[[1]])