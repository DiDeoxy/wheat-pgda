# load file paths and functions
source(file.path("src", "R", "file_paths.R"))
import::from(pgda, "calc_plot_map_stats")

# all markers
calc_plot_map_stats(
  phys_gds,
  "MAF and MR Pruned Marker Set",
  1e4,
  basename(phys_gds)
)

# ld pruned markers
calc_plot_map_stats(
  ld_gds,
  "MAF, MR, and LD Pruned Marker Set",
  1e3,
  basename(ld_gds)
)