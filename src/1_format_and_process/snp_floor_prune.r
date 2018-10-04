library(SNPRelate)

# load data
gds <- "Data\\Intermediate\\GDS\\full_phys.gds"
source("src\\R_functions\\data_loading.R")

source("src\\R_functions\\ld_floor_prune.R")

# create a list of the ld between markers on each chromosome add NA column
# and rows every million bps so that gaps appear on image
full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys.gds")
kept_ids <- by(marker_data, marker_data$chrom, function (chrom) {
    ld_mat <- snpgdsLDMat(
        full, method = "composite", snp.id = chrom$id, slide = -1
    )$LD
    # prune markers that are not in a min ld with neighbouring markers within
    # a max distance sliding window
    set.seed(1000) # so that repeated runs produce same output
    chrom$id[ld_floor_prune(chrom, ld_mat, 0.2, 10)]
})
snpgdsClose(full)
head(kept_ids[[1]])
length(unlist(kept_ids))
length(snp_id)