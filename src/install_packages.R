install.packages(
  c(
    'adegenet', 'ade', 'ape', 'BiocManager', 'circlize', 'dbscan', 'dendextend',
    'devtools', 'extrafont', 'GGally', 'ggrepel', 'igraph', 'import', 'mmod',
    'plyr', 'poppr', 'pracma', 'RColorBrewer', 'roxygen2', 'scrime', 'tidyverse'
  )
)
BiocManager::install('SNPRelate');
library(extrafont); font_import(); loadfonts()
devtools::install(file.path("src", "pgda"))
