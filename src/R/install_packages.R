install.packages(
  c(
    'adegenet', 'ade', 'ape', 'BiocManager', 'circlize', 'dbscan', 'dendextend',
    'devtools', 'emdbook', 'extrafont', 'GeneticSubsetter', 'GGally', 'ggrepel',
    'igraph', 'import', 'mmod', 'plyr', 'poppr', 'pracma', 'RColorBrewer',
    'roxygen2', 'scales', 'scrime', 'tidyverse', 'vcd'
  )
)
BiocManager::install('SNPRelate');
library(extrafont); font_import(); loadfonts()
devtools::install_github("https://github.com/DiDeoxy/PGDA")