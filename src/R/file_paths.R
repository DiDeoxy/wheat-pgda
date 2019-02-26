# raw data
raw_data <- file.path("data", "raw")
cultivars <- file.path(raw_data, "cultivars")
gene_alignments <- file.path(raw_data, "genes", "gene_alignments")
markers <- file.path(raw_data, "markers")
maps <- file.path(markers, "maps")
genotypes <- file.path(markers, "genotypes")

# formatted data
gds <- file.path("data", "gds")
ifelse(! dir.exists(gds), dir.create(gds), FALSE)
phys_gds <- file.path(gds, "maf_and_mr_pruned_phys_sample_subset.gds")
gen_gds <- file.path(gds, "maf_and_mr_pruned_gen_sample_subset.gds")
ld_gds <- file.path(gds, "ld_pruned_phys_sample_subset.gds")
geninds <- file.path("data", "geninds")
hdbscan <- file.path("data", "wheat_hdbscan.rds")
josts_d <- file.path("data", "josts_d.rds")

# results
ifelse(! dir.exists("results"), dir.create("results"), FALSE)
test_pruning <- file.path("results", "test_pruning")
ifelse(! dir.exists(test_pruning), dir.create(test_pruning), FALSE)
test_pruning_pca <- file.path(test_pruning, "pca")
ifelse(! dir.exists(test_pruning_pca), dir.create(test_pruning_pca), FALSE)
test_pruning_dend  <- file.path(test_pruning, "dend")
ifelse(! dir.exists(test_pruning_dend), dir.create(test_pruning_dend), FALSE)
josts_d_by_chrom <- file.path("results", "josts_d_by_chrom")
ifelse(! dir.exists(josts_d_by_chrom), dir.create(josts_d_by_chrom), FALSE)

# functions
functions <- file.path("src", "R", "functions")
parse_create_gds <- file.path(functions, "parse_create_gds.R")
draw_dend <- file.path(functions, "draw_dend.R")
map_stats <- file.path(functions, "map_stats.R")
locus_by_locus <- file.path(functions, "locus_by_locus.R")