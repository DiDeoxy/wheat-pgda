source("src/R_functions/funcs_calc_stats.R")

# all markers
calc_plot_map_stats(
  "phys_subset_sample",
  "Log 10 Gap Distances, Full Set",
  "LD Between Neighbouring Markers, Full Set"
)

# ld pruned markers
calc_plot_map_stats(
  "phys_subset_sample_ld_pruned",
  "Log 10 Gap Distances, Pruned Set",
  "LD Between Neighbouring Markers, Pruned Set"
)