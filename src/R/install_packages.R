# # some packages dont want to install properly on R version 3.5.3 on Fedora
# install.packages("remotes")
# library(remotes)
# install_version("fpc", "2.1-10")
# install_version("quadprog", "1.5-5")

install.packages(
  c(
    "ade4", "ape", "BiocManager", "circlize", "dbscan",
    "dendextend", "extrafont", "GeneticSubsetter", "GGally", 
    "ggrepel", "igraph", "import", "plyr", "pracma",
    "RColorBrewer", "scrime", "tidyverse", "vcd"
  )
)
library(extrafont); font_import(); loadfonts()

# have additional dependencies
install.packages(
  c(
    "adegenet", "devtools", "mmod", "poppr", "Rfast", "roxygen2"
  )
)

BiocManager::install("SNPRelate")
# devtools::install_github("https://github.com/DiDeoxy/PGDA")
devtools::install_git("../PGDA")
