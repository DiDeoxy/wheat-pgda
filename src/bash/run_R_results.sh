# create and output results

# 0
# install needed R packages
echo "install_R_packages"
Rscript src/R/install_packages.R

# 1
# calculate and plot statistics of the maf and mr pruned phys map and the ld
# pruned map
echo "calc_and_plot_map_stats"
Rscript src/R/calc_and_plot_map_stats.R

# 2
# plot a UPGMA dendrogram with hdbscan and metadata in rows of rectangels
# arranged around it
echo "plot_hdbscan_upgma_dend"
Rscript src/R/plot_hdbscan_upgma_dend.R

# 3
# calculate AMOVA phi values for various hierarchical AMOVAs
echo "calc_sample_groupings_AMOVA_phi"
Rscript src/R/calc_sample_groupings_AMOVA_phi.R

# 4
# plot the josts D values of markers between the varieties of the major
# phenptypes of the three largest clusters
echo "plot_clustered_phenos_markers_josts_ds_with_genes"
Rscript src/R/plot_clustered_phenos_markers_josts_ds_with_genes.R

# 5
# plot marker physical postions against genetic postion coloured by eh value
echo "plot_phys_vs_gen_pos_with_eh"
Rscript src/R/plot_phys_vs_gen_pos_with_eh.R

# 6
# plot ld heatmaps of each chromosome
echo "plot_LD_heatmaps"
Rscript src/R/plot_LD_heatmaps.R