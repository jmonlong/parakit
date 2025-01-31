FROM ubuntu:22.04

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

## 
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    git \
    gcc \
    make \
    bzip2 \
    tabix \
    pandoc \
    less \
    r-base r-base-dev \
    python3 \
    python3-pip \
    python3-dev \
    libxml2-dev \
    libssl-dev \
    libmariadb-dev \
    apt-transport-https \
    software-properties-common \
    dirmngr \
    gpg-agent \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    build-essential \
    pkg-config \
    libhdf5-dev \
    liblzo2-dev \
    libtokyocabinet-dev \
    libpng-dev \
    uuid-dev \
    libcurl4-openssl-dev \
    # libcurl4-gnutls-dev \
    libffi-dev \
    python3-virtualenv \
    rsync \
    python-is-python3 \
    libdeflate-dev \
    cmake \
    libjemalloc-dev \
    python3-distutils \
    pybind11-dev \
    autoconf \
    libzstd-dev \
    # libhts-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libcurl4-gnutls-dev \
    libhts-dev \
    && rm -rf /var/lib/apt/lists/*

## R
RUN R -e "install.packages(c('dplyr', 'ggplot2', 'tidyr', 'BiocManager', 'RColorBrewer', 'cowplot', 'rjson', 'rmarkdown'))" && \
    R -e "BiocManager::install(c('GenomicRanges'))"

## samtools
WORKDIR /build

RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xjvf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    ./configure && \
    make && \
    make install

## GraphAligner
ENV CONDA_DIR /opt/conda

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda install -c bioconda graphaligner

## PGGB
RUN conda install -c bioconda -c conda-forge pggb

## Minigraph-Cactus
RUN wget --no-check-certificate https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.9.3/cactus-bin-v2.9.3.tar.gz && \
    tar -xzvf cactus-bin-v2.9.3.tar.gz

# RUN cd cactus-bin-v2.9.3 && \
#     git config --global --add safe.directory /build/cactus-bin-v2.9.3 && \
#     python3 -m venv venv-cactus-v2.9.3 && \
#     printf "export PATH=$(pwd)/bin:\$PATH\nexport PYTHONPATH=$(pwd)/lib:\$PYTHONPATH\nexport LD_LIBRARY_PATH=$(pwd)/lib:\$LD_LIBRARY_PATH\n" >> venv-cactus-v2.9.3/bin/activate && \
#     chmod +x venv-cactus-v2.9.3/bin/activate && \
#     ./venv-cactus-v2.9.3/bin/activate && \
#     python3 -m pip install -U setuptools pip wheel && \
#     python3 -m pip install -U . && \
#     python3 -m pip install -U -r ./toil-requirement.txt

ENV PATH=/build/cactus-bin-v2.9.3/bin:$PATH
ENV PYTHONPATH=/build/cactus-bin-v2.9.3/lib:$PYTHONPATH
ENV LD_LIBRARY_PATH=/build/cactus-bin-v2.9.3/lib:$LD_LIBRARY_PATH

RUN cd cactus-bin-v2.9.3 && \
    git config --global --add safe.directory /build/cactus-bin-v2.9.3 && \
    python3 -m pip install -U setuptools pip wheel && \
    python3 -m pip install -U . && \
    python3 -m pip install -U -r ./toil-requirement.txt

## vg
ENV PATH=$PATH:/bin

WORKDIR /bin
RUN wget --no-check-certificate https://github.com/vgteam/vg/releases/download/v1.63.1/vg && \
    chmod +x vg

## AGC https://github.com/refresh-bio/agc/releases
RUN wget --no-check-certificate https://github.com/refresh-bio/agc/releases/download/v3.2.1/agc-3.2_x64_linux.tar.gz && \
    tar -xzvf agc-3.2_x64_linux.tar.gz && \
    mv agc-3.2_x64_linux/agc . && \
    chmod +x agc && \
    rm -r agc-3.2_x64_linux agc-3.2_x64_linux.tar.gz

RUN pip3 install snakemake

## parakit
COPY . /parakit

WORKDIR /parakit

RUN python3 -m pip install -e .

WORKDIR /snakemake_cache

RUN chmod a+xrw /snakemake_cache

ENV XDG_CACHE_HOME=/snakemake_cache

WORKDIR /home

