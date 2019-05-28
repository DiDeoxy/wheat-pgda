# format and process data R

# 0
# install needed R packages
echo "install_R_packages"
Rscript src/R/install_packages.R

# 1
# take the two genetic maps we have and compare them to identify markers
# assigned to the same linkage groups
echo "create_filtered_map"
Rscript src/R/create_filtered_map.R

# 2
# create a tibble using the consensus markers from step 1 with the physical
# positions from the gmap alignment of the probes supplied by the IWGSC
# with the genotype information for those markers supplied by Dr. Pozniak
echo "create_phys_map_genotypes_tibble"
Rscript src/R/create_phys_map_genotypes_tibble.R

# 3
# create gds format files from the table created in step 2 one with gmap
# postions, the other with pozniak genetic map positions
echo "create_gds"
Rscript src/R/create_gds.R

# 4
# create a series of gds files, first pruning sets of samples with >99% IBS
# then removing markers with a MAF < 0.05 or a MR > 0.10, finally LD prune
echo "filter_maf_mr_sample_and_ld_prune"
Rscript src/R/filter_maf_mr_sample_and_ld_prune.R

# 5
# calculate HDBSCAN clusters
echo "hdbscan_cluster"
Rscript src/R/hdbscan_cluster.R

# 6
# create geninds for calculating Jost's Ds and AMOVA phis between groups
echo "create_geninds"
Rscript src/R/create_geninds.R

# 7
# calculate Jost's D values
echo "calc_josts_d"
Rscript src/R/calc_josts_d.R