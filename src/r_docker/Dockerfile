################# BASE IMAGE ######################
FROM rocker/r-base

################## MAINTAINER #####################
LABEL maintainer.name="William Hargreaves"
LABEL maintainer.email="whargrea@uoguelph.ca"

################# INSTALLATION ####################
RUN BUILDDEPS="\
  libssl-dev \
  libxml2-dev \
  libcurl4-openssl-dev \
  libpq-dev \
  libmysql++-dev \
  libudunits2-dev \
  libgdal-dev \
  default-jdk \
  r-cran-rjava" \
  && apt-get update \
  && apt-get -y --no-install-recommends install $BUILDDEPS \
  && install2.r --error --deps TRUE \
  tidyverse \
  dplyr \
  plyr \
  rJava \
  GGally \
  ggrepel \
  ape \
  circlize \
  dendextend \
  extrafont \
  scrime \
  dbscan \
  adegenet \
  stringr \
  poppr \
  RColorBrewer \
  pracma \
  BiocManager \
  && Rscript -e "BiocManager::install('SNPRelate')" \
  && Rscript -e "library(extrafont); font_import(); loadfonts()" \
  && apt-get remove --purge -y $BUILDDEPS \
  && apt-get autoremove -y \
  && apt-get autoclean -y \
  && rm -rf /var/lib/apt/lists/*