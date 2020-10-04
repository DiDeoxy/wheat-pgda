# load file paths and functions
source("wheat-pgda/src/R/file_paths.R")
import::from(pgda, "calc_plot_map_stats")

# all markers
calc_plot_map_stats(
  phys_gds,
  gen_gds,
  "full_set",
  "Full Marker Set Histograms",
  1e4
)

# all markers
calc_plot_map_stats(
  ld_phys_gds,
  ld_gen_gds,
  "ld_pruned_set",
  "LD Pruned Marker Set Histograms",
  1e3
)