# load functions
source(map_stats)

# all markers
calc_plot_map_stats(
  phys_gds,
  "Log 10 Gap Distances, Full Set",
  "LD Between Neighbouring Markers, Full Set",
  1500
)

# ld pruned markers
calc_plot_map_stats(
  ld_gds,
  "Log 10 Gap Distances, Pruned Set",
  "LD Between Neighbouring Markers, Pruned Set",
  500
)