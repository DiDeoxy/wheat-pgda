library(SNPRelate)

# apply the best LD pruning to the dataset
gds <- "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
source("src\\R_functions\\data_loading.R")

# performs LD pruning to produce a set of SNPs that maximally represent the
# diversity of the genome with as little redundant info as possible used for
# estimation of relationships genome wide (i.e. clustering)
full <- snpgdsOpen("Data\\Intermediate\\GDS\\full_phys_subset_sample.gds")
set.seed(1000)
snpset_ids_list <- snpgdsLDpruning(
    full,
    autosome.only = F, missing.rate = 0.1, maf = 0.05,
    ld.threshold = 0.6, slide.max.bp = 1e7
)
snpgdsClose(full)
snp_indices <- match(unlist(snpset_ids_list), snp_id)

# create a subset gds object containg only the SNPs remaing after ld pruning
snpgdsCreateGeno("Data\\Intermediate\\GDS\\full_phys_subset_sample_pruned.gds",
    genmat = genotypes[snp_indices, ],
    sample.id = sample_id,
    snp.id = snp_id[snp_indices],
    snp.chromosome = snp_chrom[snp_indices],
    snp.position = snp_pos[snp_indices],
    other.vars = list(samp_annot = samp_annot),
    snpfirstdim = T
)