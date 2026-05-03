FROM ubuntu:24.04

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    git \
    wget \
    rsync \
    cmake \
    autoconf \
    python3 \
    python3-dev \
    python3-pip \
    python3.12-venv \
    python3-setuptools \
    r-base r-base-dev \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

##
## Start with R (long to install)
##

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

##
## Cactus, also quite long to build
##

# RUN git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive

RUN git clone https://github.com/ComparativeGenomicsToolkit/cactus.git && \
    cd cactus && git checkout cd625082f6038d64b4a1e4762cd2fde50cc2f657 && git submodule update --init --recursive

# Create and activate the virtual environment
RUN python3 -m venv /opt/venv
ENV PATH "/opt/venv/bin:$PATH"

# # Upgrade pip and setuptools
# RUN pip install --upgrade pip setuptools wheel

WORKDIR /build/cactus

# ENV PATH=/build/cactus/bin:$PATH
# ENV PYTHONPATH=/build/cactus/lib:$PYTHONPATH
# ENV LD_LIBRARY_PATH=/build/cactus/lib:$LD_LIBRARY_PATH

RUN python3 -m pip install -U setuptools pip wheel && \
    python3 -m pip install -U . && \
    python3 -m pip install -U -r ./toil-requirement.txt

RUN apt-get update && apt-get install -y \
    build-essential git python3 python3-dev python3-pip zlib1g-dev wget libbz2-dev pkg-config libhdf5-dev liblzo2-dev libtokyocabinet-dev wget liblzma-dev libxml2-dev libssl-dev libpng-dev uuid-dev libcurl4-gnutls-dev libffi-dev python3-virtualenv rsync python-is-python3 libdeflate-dev cmake libjemalloc-dev pybind11-dev autoconf libzstd-dev liblz4-dev libhts-dev

RUN make

RUN build-tools/downloadPangenomeTools

## test with
## cd /build/cactus
## cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal

##
## Tools in Conda: GraphAligner, PGGB
##

ENV CONDA_DIR /opt/conda

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

## GraphAligner
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda install -c bioconda graphaligner

## PGGB
RUN conda install -c bioconda -c conda-forge pggb


##
## Tools with binaries
##

ENV PATH=$PATH:/bin
WORKDIR /bin

## vg
RUN wget --no-check-certificate https://github.com/vgteam/vg/releases/download/v1.63.1/vg && \
    chmod +x vg

## AGC https://github.com/refresh-bio/agc/releases
RUN wget --no-check-certificate https://github.com/refresh-bio/agc/releases/download/v3.2.1/agc-3.2_x64_linux.tar.gz && \
    tar -xzvf agc-3.2_x64_linux.tar.gz && \
    mv agc-3.2_x64_linux/agc . && \
    chmod +x agc && \
    rm -r agc-3.2_x64_linux agc-3.2_x64_linux.tar.gz

## minimap2
RUN wget --no-check-certificate https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30_x64-linux.tar.bz2 && \
    tar -xjvf minimap2-2.30_x64-linux.tar.bz2 && \
    mv minimap2-2.30_x64-linux/minimap2 . && \
    chmod +x minimap2 && \
    rm -r minimap2-2.30_x64-linux.tar.bz2 minimap2-2.30_x64-linux

##
## Parakit
##

COPY . /parakit

WORKDIR /parakit

RUN python3 -m pip install -e .


## also install snakemake
RUN pip3 install snakemake

WORKDIR /snakemake_cache

RUN chmod a+xrw /snakemake_cache

ENV XDG_CACHE_HOME=/snakemake_cache

## end at home
WORKDIR /home
