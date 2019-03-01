# functions
functions <- file.path("src", "functions")
snpgds_parse_subset <- file.path(functions, "snpgds_parse_subset.R")
draw_dend <- file.path(functions, "draw_dend.R")
map_stats <- file.path(functions, "map_stats.R")
locus_by_locus <- file.path(functions, "locus_by_locus.R")

################################################################################

# raw data
raw_data <- file.path("data", "raw")
cultivars <- file.path(raw_data, "cultivars")
markers <- file.path(raw_data, "markers")
gene_alignments <- file.path(raw_data, "gene_alignments")

################################################################################
# intermediate data
intermediate <- file.path("data", "intermediate")
ifelse(! dir.exists(intermediate), dir.create(intermediate), FALSE)

inter_markers <- file.path(intermediate, "markers")
ifelse(! dir.exists(inter_markers), dir.create(inter_markers), FALSE)

gds <- file.path(intermediate, "gds")
ifelse(! dir.exists(gds), dir.create(gds), FALSE)

phys_gds <- file.path(gds, "maf_and_mr_pruned_phys_sample_subset.gds")
gen_gds <- file.path(gds, "maf_and_mr_pruned_gen_sample_subset.gds")
ld_gds <- file.path(gds, "ld_pruned_phys_sample_subset.gds")

geninds <- file.path(intermediate, "geninds")
ifelse(! dir.exists(geninds), dir.create(geninds), FALSE)

hdbscan <- file.path(intermediate, "wheat_hdbscan.rds")
josts_d <- file.path(intermediate, "josts_d.rds")

################################################################################
# results
ifelse(! dir.exists("results"), dir.create("results"), FALSE)
map_stats_and_plots <- file.path("results", "map_stats_and_plots")
ifelse(
  ! dir.exists(map_stats_and_plots), dir.create(map_stats_and_plots), FALSE
)
josts_d_by_chrom <- file.path("results", "josts_d_by_chrom")
ifelse(! dir.exists(josts_d_by_chrom), dir.create(josts_d_by_chrom), FALSE)