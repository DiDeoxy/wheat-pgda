################################################################################
# install packages
source("https://bioconductor.org/biocLite.R")
biocLite("SNPRelate")

install.packages(
  c(
    "tidyverse", "plyr", "GGally", "ggrepel", "ape", "adegenet", "mmod",
    "poppr", "dbscan", "scrime", "circlize", "dendextend", "RColorBrewer",
    "extrafont"
  )
)

library(extrafont)
font_import()
loadfonts()

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

# load functions
source(gds_parse_create)

################################################################################
# format and process data

# 1
source(file.path(base "generate_filtered_maps.R"))

# 2
library(ape)
source(file.path(base, "create_phys_map_genotypes_tibble.R"))

# 3
source(file.path(base, "create_gds.R"))

# 4
source(file.path(base, "remove_NILs_and_prune_markers.R"))

# 5
library(circlize)
source(draw_dend)
library(dendextend)
library(GGally)
library(extrafont)
source(file.path(base, "test_ld_pruning_parameters.R"))

# 6
source(file.path(base, "ld_prune.R"))

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
source(calc_map_stats)
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
source(file.path(base, "plot_clustered_pheno_all_genes.R"))

# 5
source(file.path(base, "plot_phys_vs_gen_pos.R"))

# 6
source(file.path(base, "plot_LD_heatmaps.R"))
