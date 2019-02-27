# install needed packages
# Rscript src/R/install.R

################################################################################
# format and process data

# 1
# take the two genetic maps we have and compare them to identify markers
# assigned to the same linkage groups
echo "generate_filtered_maps"
Rscript src/generate_filtered_maps.R

# 2
# create a tibble using the consensus markers from step 1 with the physical
# positions from the gmap alignment of the probes supplied by the IWGSC
# with the genotype information for those markers supplied by Dr. Pozniak
echo "create_phys_map_genotypes_tibble"
Rscript src/create_phys_map_genotypes_tibble.R

# 3
# create gds format files from the table created in step 2 one with gmap
# postions, the other with pozniak genetic map positions
echo "create_gds"
Rscript src/create_gds.R

# 4
# create a series of gds files, first pruning sets of samples with >99% IBS
# then removing markers with a MAF < 0.05 or a MR > 0.10, finally LD prune
echo "sample_maf_mr_and_ld_prune_data"
Rscript src/sample_maf_mr_and_ld_prune_data.R

# 5
# calculate HDBSCAN clusters
echo "hdbscan_cluster"
Rscript src/hdbscan_cluster.R

# 6
# create geninds for calculating Jost's Ds and AMOVA phis between groups
echo "create_geninds"
Rscript src/create_geninds.R

# 7
# calculate Jost's D values
echo "calc_josts_d"
Rscript src/calc_josts_d.R

################################################################################
# create and output results

# 1
# calculate and plot statistics of the maf and mr pruned phys map and the ld
# pruned map
echo "calc_and_plot_map_stats"
Rscript src/calc_and_plot_map_stats.R