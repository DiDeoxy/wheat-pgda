# install needed packages
Rscript src/install.R

################################################################################
# format and process data

# 1
# take the two genetic maps we have and compare them to identify markers
# assigned to the same linkage groups
echo "create_filtered_map"
Rscript src/create_filtered_map.R

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

# 2
# plot a UPGMA dendrogram with hdbscan and metadata in rows of rectangels
# arranged around it
echo "plot_hdbscan_upgma_dend"
Rscript src/plot_hdbscan_upgma_dend.R

# 3
# calculate AMOVA phi values for various hierarchical AMOVAs
echo "calc_sample_groupings_AMOVA_phi"
Rscript src/calc_sample_groupings_AMOVA_phi.R

#4
# plot the josts D values of markers between the varieties of the major
# phenptypes of the three largest clusters
echo "plot_clustered_phenos_markers_josts_ds_with_genes"
Rscript src/plot_clustered_phenos_markers_josts_ds_with_genes.R

# 5
# plot marker physical postions against genetic postion coloured by eh value
echo "plot_phys_vs_gen_pos_with_eh"
Rscript src/plot_phys_vs_gen_pos_with_eh.R

# 6
# plot ld heatmaps of each chromosome
echo "plot_LD_heatmaps"
Rscript src/plot_LD_heatmaps.R