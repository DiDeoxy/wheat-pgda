# # some packages dont want to install properly on R version 3.5.3 on Fedora
# install.packages("remotes")
# library(remotes)
# install_version("fpc", "2.1-10")
# install_version("quadprog", "1.5-5")

install.packages(
  c(
    "adegenet", "ade4", "ape", "BiocManager", "circlize", "dbscan",
    "dendextend", "devtools", "extrafont", "GeneticSubsetter", "GGally", 
    "ggrepel", "igraph", "import", "mmod", "plyr", "poppr", "pracma",
    "RColorBrewer", "roxygen2", "scrime", "tidyverse", "vcd"
  )
)
library(extrafont); font_import(); loadfonts()


BiocManager::install("SNPRelate")
devtools::install_github("https://github.com/DiDeoxy/PGDA")