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

blast <- file.path(intermediate, "blast")

hdbscan_rds <- file.path(intermediate, "wheat_hdbscan.rds")
josts_d_rds <- file.path(intermediate, "josts_d.rds")

################################################################################
# results
ifelse(! dir.exists("results"), dir.create("results"), FALSE)

map_stats_and_plots <- file.path("results", "map_stats_and_plots")
ifelse(
  ! dir.exists(map_stats_and_plots), dir.create(map_stats_and_plots), FALSE
)

allele_richness <- file.path("results", "allele_richness")
ifelse(! dir.exists(allele_richness), dir.create(allele_richness), FALSE)

zoomed_marker_plots <- file.path("results", "zoomed_marker_plots")
ifelse(
  ! dir.exists(zoomed_marker_plots), dir.create(zoomed_marker_plots), FALSE
)

all_data_chroms <- file.path("results", "all_data_chroms")
ifelse(
  ! dir.exists(all_data_chroms), dir.create(all_data_chroms), FALSE
)