source("src\\R_functions\\funcs_calc_map_stats.R")

# all markers
calc_map_stats_plot(
  "phys_subset_sample",
  "Log 10 LD Gap Distances All Markers"
)

# ld pruned markers
calc_map_stats_plot(
  "phys_subset_sample_ld_pruned",
  "Log 10 LD Gap Distances Ld Pruned Markers"
)