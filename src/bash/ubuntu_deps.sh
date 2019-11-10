# # install latest R (need to do manually)
# sudo apt install -y apt-transport-https software-properties-common
# sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
# sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
# sudo apt-get update
# sudo apt-get install -y r-base r-base-dev

sudo apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libudunits2-dev \
    libgdal-dev \
    libgsl-dev \
    libxml2-dev