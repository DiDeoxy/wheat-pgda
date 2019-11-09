# debian buster
FROM debian:10.1

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# update apt-get and install needed packages
RUN apt-get update --fix-missing \
    && apt-get install -y --no-install-recommends \
        sudo \
        wget \
        bzip2 \
        ca-certificates \
        curl \
        bash-completion \
        less \
        ssh \
        git \
        pandoc \
        libcairo2-dev \
        libopenblas-dev \
        nano \
        # install git iproute2, process tools, lsb-release (common in install
        # instructions for CLIs) and libzip for R Tools extension
        iproute2 \
        procps \
        lsb-release \
        libzip-dev \
    # clean up
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# install .Net runtime for R Tools extension
# this doesn't work when attached to above (something to do with wget)
RUN wget --quiet https://packages.microsoft.com/config/debian/10/packages-microsoft-prod.deb \
    && dpkg -i packages-microsoft-prod.deb \
    && apt-get install apt-transport-https \
    && apt-get update \
    && apt-get -y install dotnet-runtime-2.2 \
    # clean up again
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

ARG USERNAME=maxh
# On Linux, replace with your actual UID, GID if not the default 1000
ARG USER_UID=1001
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && mkdir -p /home/$USERNAME/.vscode-server /home/$USERNAME/.vscode-server-insiders \
    && chown ${USER_UID}:${USER_GID} /home/$USERNAME/.vscode-server* \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Set the default user
USER $USERNAME
WORKDIR /home/$USERNAME

# setup ssh
RUN echo 'eval $(ssh-agent -s)' >> /home/$USERNAME/.bashrc \
    && echo 'ssh-add ~/.ssh/id_rsa' >> /home/$USERNAME/.bashrc

# install conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O /home/maxh/miniconda.sh \
    && /bin/bash /home/$USERNAME/miniconda.sh -b -p /home/$USERNAME/conda \
    && rm /home/$USERNAME/miniconda.sh \
    && /home/$USERNAME/conda/bin/conda clean -tipsy \
    && sudo ln -s /home/$USERNAME/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
    && echo ". /home/$USERNAME/conda/etc/profile.d/conda.sh" >> /home/$USERNAME/.bashrc \
    && echo "conda activate base" >> /home/$USERNAME/.bashrc

# add conda to path
ENV PATH /home/$USERNAME/conda/bin:$PATH

# install R through conda
RUN conda update conda \
    && conda install -c conda-forge r

# set default R mirror
RUN echo '# Default repo \n\
local({ \n\
    r <- getOption("repos") \n\
    r["CRAN"] <- "https://cloud.r-project.org" \n\
    options(repos=r) \n\
})' >> /home/maxh/.Rprofile

# install R packages available on conda-forge
RUN conda install -c conda-forge \
        r-ade4 \
        r-ape \
        r-bh \
        r-biocmanager \
        r-circlize \
        r-dendextend \
        r-devtools \
        r-extrafont \
        r-ggally \
        r-ggrepel \
        r-ggpubr \
        r-igraph \
        r-import \
        r-mmod \
        r-phangorn \
        r-plyr \
        r-pracma \
        r-quadprog \
        r-rcolorbrewer \
        r-rfast \
        r-roxygen2 \
        r-scrime \
        r-tidyverse \
        r-vcd \
    && conda install -c bioconda \
        r-adegenet \
        bioconductor-snprelate \
    && conda install -c anaconda make gcc_linux-64 

RUN sudo ln -s /bin/tar /bin/gtar \
    && Rscript -e 'install.packages(c("stringi", "dbscan", "GeneticSubsetter", "poppr"))' \
    && Rscript -e 'devtools::install_github("DiDeoxy/pgda", dep = FALSE)'

# install powerline
RUN pip install powerline-status \
    && echo "powerline-daemon -q" >> /home/$USERNAME/.bashrc \
    && echo ". /home/maxh/conda/lib/python3.7/site-packages/powerline/bindings/bash/powerline.sh"  >> /home/$USERNAME/.bashrc

# install tini?
ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN sudo chmod +x /usr/bin/tini

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]