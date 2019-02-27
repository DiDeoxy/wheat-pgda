################################################################################
# load colours
source(file.path("src", "colour_sets.R"))

# load common functions
source(parse_create_gds)

################################################################################
# format and process data

# load the needed libraries and file paths
suppressPackageStartupMessages(library(tidyverse))
library(SNPRelate)
source(file.path("src", "file_paths.R"))
source(parse_create_gds)

################################################################################
# output results

# 1


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
