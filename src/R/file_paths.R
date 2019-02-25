# raw data
markers <- file.path("data", "markers")
maps <- file.path(markers, "maps")
genotypes <- file.path(markers, "genotypes")
cultivars <- file.path("data", "cultivars")
genes <- file.path("data", "genes")
selected_alignments <- file.path(genes, "gene_alignments")

# formatted data
gds <- file.path("data", "gds")
ifelse(! dir.exists(gds), dir.create(gds), FALSE)
phys_gds <- file.path(gds "maf_and_mr_pruned_phys_sample_subset.gds")
gen_gds <- file.path(gds "maf_and_mr_pruned_gen_sample_subset.gds")
ld_gds <- file.path(gds "ld_pruned_phys_sample_subset.gds")
geninds <- file.path("data", "geninds")
josts_d <- file.path("data", "josts_d.rds")

# results
ifelse(! dir.exists("results"), dir.create("results"), FALSE)
test_pruning <- file.path("results", "test_pruning")
test_pruning_pca <- file.path(test_pruning, "pca")
ifelse(! dir.exists(test_pruning_pca), dir.create(test_pruning_pca), FALSE)
test_pruning_dend  <- file.path(test_pruning, "dend")
ifelse(! dir.exists(test_pruning_dend), dir.create(test_pruning_dend), FALSE)
hdbscan <- file.path("results", "wheat_hdbscan.rds")
adegenet <- file.path("results", "adegenet")
ifelse(! dir.exists(adegenet), dir.create(adegenet), FALSE)
dend <- file.path("results", "dend.png")

# functions
functions <- file.path("src", "R", "functions")
gds_parse_create <- file.path(functions "gds_parse_create.R")
draw_dend <- file.path(functions, "draw_dend.R")
calc_plot_map_stats <- file.path(functions, "calc_plot_map_stats.R")
locus_by_locus <- file.path(functions, "locus_by_locus.R")