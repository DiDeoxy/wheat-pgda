install.packages("extrafont")
library(extrafont)
font_import()
loadfonts()

source("https://bioconductor.org/biocLite.R")
biocLite("SNPRelate")

install.packages(
  c(
    "tidyverse", "plyr", "GGally", "ggrepel", "ape", "adegenet", "mmod",
    "poppr", "dbscan", "scrime", "circlize", "dendextend", "RColorBrewer"
  )
)