# raw data
raw_data <- file.path("data", "raw")
cultivars <- file.path(raw_data, "cultivars")
markers <- file.path(raw_data, "markers")

################################################################################
# intermediate data
intermediate <- file.path("data", "intermediate")
ifelse(! dir.exists(intermediate), dir.create(intermediate), FALSE)

inter_markers <- file.path(intermediate, "markers")
ifelse(! dir.exists(inter_markers), dir.create(inter_markers), FALSE)

gds <- file.path(intermediate, "gds")
ifelse(! dir.exists(gds), dir.create(gds), FALSE)

phys_gds <- file.path(gds, "sample_pruned_phys.gds")
gen_gds <- file.path(gds, "sample_pruned_gen.gds")
ld_phys_gds <- file.path(gds, "maf_mr_filtered_sample_ld_pruned_phys.gds")
ld_gen_gds <- file.path(gds, "maf_mr_filtered_sample_ld_pruned_gen.gds")

geninds <- file.path(intermediate, "geninds")
ifelse(! dir.exists(geninds), dir.create(geninds), FALSE)

hdbscan <- file.path(intermediate, "wheat_hdbscan.rds")
josts_d <- file.path(intermediate, "josts_d.rds")

blast <- file.path(intermediate, "blast")

################################################################################
# results
ifelse(! dir.exists("results"), dir.create("results"), FALSE)

map_stats_and_plots <- file.path("results", "map_stats_and_plots")
ifelse(
  ! dir.exists(map_stats_and_plots), dir.create(map_stats_and_plots), FALSE
)

josts_d_by_chrom <- file.path("results", "josts_d_by_chrom")
ifelse(! dir.exists(josts_d_by_chrom), dir.create(josts_d_by_chrom), FALSE)

allele_richness <- file.path("results", "allele_richness")
ifelse(! dir.exists(allele_richness), dir.create(allele_richness), FALSE)