library(SNPRelate)

ld_prune_densest_first <- function(params, subset) {
  source(str_c("src\\R_functions\\funcs_ld_prune_", subset, ".r"))
  subset_sample <- snpgdsOpen(
    "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
  )
  results <- unlist(
    ld_prune(
      subset_sample, params$maf, params$window_size, params$ld_floor,
      params$ld_ceiling
    )
  )
  snpgdsClose(subset_sample)
  results
}

ld_prune_max_dist_max_snps <- function(params, subset, snp_id) {
  source(str_c("src\\R_functions\\funcs_ld_prune_", subset, ".r"))
  subset_sample <- snpgdsOpen(
    "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
  )
  kept_id <- unlist(
    ld_prune(
      subset_sample, params$maf, params$max_snps, params$max_dist,
      params$ld_floor, params$ld_ceiling
    )
  )
  snpgdsClose(subset_sample)
  match(kept_id, snp_id)
}

ld_prune_max_dist <- function(params, subset, snp_id) {
  source(str_c("src\\R_functions\\funcs_ld_prune_", subset, ".r"))
  subset_sample <- snpgdsOpen(
    "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
  )
  kept_id <- unlist(
    ld_prune(
      subset_sample, params$maf, params$max_dist, params$ld_floor,
      params$ld_ceiling
    )
  )
  snpgdsClose(subset_sample)
  match(kept_id, snp_id)
}

ld_prune_max_snps <- function(params, subset, snp_id) {
  source(str_c("src\\R_functions\\funcs_ld_prune_", subset, ".r"))
  subset_sample <- snpgdsOpen(
    "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
  )
  kept_id <- unlist(
    ld_prune(
      subset_sample, params$maf, params$max_snps, params$ld_floor,
      params$ld_ceiling
    )
  )
  snpgdsClose(subset_sample)
  match(kept_id, snp_id)
}

ld_prune_scaled_snps <- function(params, subset, snp_id) {
  source(str_c("src\\R_functions\\funcs_ld_prune_", subset, ".r"))
  subset_sample <- snpgdsOpen(
    "Data\\Intermediate\\GDS\\full_phys_subset_sample.gds"
  )
  kept_id <- unlist(
    ld_prune(
      subset_sample, params$maf, params$min_snps, params$max_snps,
      params$ld_floor, params$ld_ceiling
    )
  )
  snpgdsClose(subset_sample)
  match(kept_id, snp_id)
}