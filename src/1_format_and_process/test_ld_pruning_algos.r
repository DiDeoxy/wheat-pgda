library(SNPRelate)
library(tidyverse)

source("src\\R_functions\\funcs_gds_parse_create.R")
source("src\\R_functions\\funcs_calc_stats.R")
source("src\\R_functions\\funcs_plot_eh.R")
source("src\\R_functions\\funcs_call_prune.R")

params1 <- list(
  maf = 0.05, min_snps = 3, max_snps = 21, max_dist = 1e7,
  ld_floor = 0.2, ld_ceiling = 0.5, window_size = 20
)
# params2 <- list(
#   maf = 0, min_snps = 3, max_snps = 21, max_dist = 1e7, ld_floor = 0.2,
#   ld_ceiling = 1
# )

subsets <- c(
  "max_dist_max_snps", "max_dist", "max_snps", "scaled_snps"
)

results <- ld_prune_densest_first(params1, "densest_first")

# subset_lengths <- ""
# for (params in list(params1, params2)[1]) {
#   indices <- list()
#   for (i in 1:length(subsets)) {
#     wheat <- load_data("full_phys_subset_sample", mb = FALSE)
#     if (subsets[i] == "max_dist_max_snps")  {
#       # print(i)
#       indices[[i]] <- ld_prune_max_dist_max_snps(
#         params, subsets[i], wheat$snp_data$id)
#     } else if (subsets[i] == "max_dist") {
#       # print(i)
#       indices[[i]] <- ld_prune_max_dist(
#         params, subsets[i], wheat$snp_data$id)
#     } else if (subsets[i] == "max_snps")  {
#       # print(i)
#       indices[[i]] <- ld_prune_max_snps(
#         params, subsets[i], wheat$snp_data$id)
#     } else if (subsets[i] == "scaled_snps") {
#       # print(i)
#       indices[[i]] <- ld_prune_scaled_snps(
#         params, subsets[i], wheat$snp_data$id)
#     }
#     snpgds_create_snp_subset(wheat, indices[[i]], subsets[i])

#     wheat <- load_data(subsets[i], mb = FALSE)
#     calc_map_stats_plot(wheat, subsets[i])
#     plot_eh(wheat, subsets[i])

#     subset_lengths <- str_c(
#       subset_lengths, subsets[i], ": ", length(indices[[i]]), "\n"
#     )
#   }
# }

# cat(
#   "maf of ", params$maf, " and  ld_ceiling, of ",
#   params$ld_ceiling, "\n", subset_lengths
# )

# start with densets regions find center of mass find reasonable window size
# from outside that central point, find the clique in the region and remove all
# markers encompassed by the region of the clique

# desnity is defined by the average distance to the nearest 20 snps

# find largest cliques in the region iteratively and chenck for correlation,
# if they are uncorrelated be suspicious if correlated take the largest