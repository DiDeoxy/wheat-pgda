install.packages(
  c(
    'adegenet', 'ape', 'BiocManager', 'circlize', 'dendextend', 'dbscan', 
    'extrafont', 'GGally', 'ggrepel', 'igraph', 'import', 'mmod', 'plyr', 'poppr',
    'pracma', 'RColorBrewer', 'scrime', 'tidyverse' 
  )
)
BiocManager::install('SNPRelate');
library(extrafont); font_import(); loadfonts()
