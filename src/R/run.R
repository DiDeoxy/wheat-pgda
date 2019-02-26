################################################################################
# install packages
# source("https://bioconductor.org/biocLite.R")
# biocLite("SNPRelate")

# install.packages(
#   c(
#     "tidyverse", "plyr", "GGally", "ggrepel", "ape", "adegenet", "mmod",
#     "poppr", "dbscan", "scrime", "circlize", "dendextend", "RColorBrewer",
#     "extrafont", "pracma"
#   )
# )

# library(extrafont)
# font_import()
# loadfonts()

################################################################################
# load libraries
library(tidyverse)
library(SNPRelate)

# base file path for scripts
base <- file.path("src", "R")

# load data and results paths
source(file.path(base, "file_paths.R"))

# load colours
source(file.path(base, "colour_sets.R"))

# load common functions
source(parse_create_gds)

################################################################################
# format and process data

# 1
source(file.path(base, "generate_filtered_maps.R"))

# 2
library(ape)
source(file.path(base, "create_phys_map_genotypes_tibble.R"))

# 3
source(file.path(base, "create_gds.R"))

# 4
source(file.path(base, "sample_maf_mr_and_ld_prune_data.R"))

# 7
library(scrime)
library(dbscan)
source(file.path(base, "hdbscan_cluster.R"))

# 8
library(adegenet)
library(plyr)
source(file.path(base, "create_geninds.R"))

# 9
library(mmod)
source(file.path(base, "calc_josts_d.R"))

################################################################################
# output results

# 1
library(ggplot2)
library(GGally)
library(extrafont)
library(RColorBrewer)
library(pracma)
source(map_stats)
source(file.path(base, "calc_and_plot_map_stats.R"))

# 2
library(circlize)
source(draw_dend)
library(dendextend)
library(extrafont)
library(ape)
source(file.path(base, "plot_hdbscan_upgma_dend.R"))

# 3
library(poppr)
library(ape)
library(stringr)
source(file.path(base, "calc_sample_groupings_AMOVA_phi.R"))

# 4
library(plyr)
library(GGally)
library(ggrepel)
library(extrafont)
source(locus_by_locus)
source(file.path(base, "plot_clustered_phenos_markers_josts_ds_with_genes.R"))

# these two arent outputting for some reason
# 5
library(plyr)
library(GGally)
library(extrafont)
source(locus_by_locus)
source(file.path(base, "plot_phys_vs_gen_pos_with_eh.R"))

# 6
library(RColorBrewer)
library(extrafont)
source(file.path(base, "plot_LD_heatmaps.R"))
