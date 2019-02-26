################# BASE IMAGE ######################
FROM debian:stretch-slim

################## MAINTAINER #####################
LABEL license="GPL-2.0" \
  maintainer.name="William Hargreaves" \
  maintainer.email="whargrea@uoguelph.ca"

################# INSTALLATION ####################
# setup environment
ENV BUILD_DATE=2019-02-26 \
  R_VERSION=3.5.2 \
  LANGUAGE=en_US.UTF-8 \
  LC_ALL=en_US.UTF-8 \
  LANG=en_US.UTF-8 \
  TERM=xterm

# set locale
RUN apt-get update \
  && apt-get install -y --no-install-recommends locales \
  && echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

# install build dependencies
RUN BUILDDEPS="\
  curl \
  libcurl3 \
  # zlib1g \
  make " \
  && apt-get install -y --no-install-recommends $BUILDDEPS

# Download, make, and install R
RUN cd tmp/ \
  && curl -O http://cran.utstat.utoronto.ca/src/base/R-3/R-${R_VERSION}.tar.gz \
  tar -xf R-${R_VERSION}.tar.gz \
  && cd R-${R_VERSION} \
  ## Set compiler flags
  && R_PAPERSIZE=letter \
  R_BATCHSAVE="--no-save --no-restore" \
  R_BROWSER=xdg-open \
  PAGER=/usr/bin/pager \
  PERL=/usr/bin/perl \
  R_UNZIPCMD=/usr/bin/unzip \
  R_ZIPCMD=/usr/bin/zip \
  R_PRINTCMD=/usr/bin/lpr \
  LIBnn=lib \
  AWK=/usr/bin/awk \
  CFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
  CXXFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
  ## Configure options
  ./configure --enable-R-shlib \
  --enable-memory-profiling \
  --with-readline \
  --with-blas \
  --with-tcltk \
  --disable-nls \
  --with-recommended-packages \
  ## Build and install
  && make \
  && make install