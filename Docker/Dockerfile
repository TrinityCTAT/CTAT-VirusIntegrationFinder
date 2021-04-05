FROM ubuntu:18.04
MAINTAINER bhaas@broadinstitute.org
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    g++ \
    gcc \
    git \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libjpeg-dev \
    liblzma-dev \
    libpng-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    make \
    pbzip2 \
    pigz \
    python3 \
    python3-pip \
    r-base \
    unzip \
    wget \
    x11-common \
    zlib1g \
    zlib1g-dev \
    zlibc

# make python3 be the default python
RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN pip3 install igv-reports==1.0.2 pandas==1.1.5

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("tidyverse")'
RUN R -e 'BiocManager::install("argparse")'

## Samtools
WORKDIR $SRC
ENV SAMTOOLS_VERSION=1.9
RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    cd $SRC && \
    wget $SAMTOOLS_URL && \
    tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    rm $SRC/samtools-${SAMTOOLS_VERSION}.tar.bz2 && rm -r $SRC/samtools-${SAMTOOLS_VERSION}

## Picard tools
WORKDIR $SRC
RUN wget https://github.com/broadinstitute/picard/releases/download/2.25.1/picard.jar
ENV PICARD_HOME $SRC/picard.jar

## STAR Aligner
WORKDIR $SRC
ENV STAR_VERSION=2.7.8a
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar xf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin && \
    rm $SRC/${STAR_VERSION}.tar.gz && rm -r $SRC/STAR-${STAR_VERSION}

## Bedtools
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
   tar xf bedtools-2.30.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/ && \
   rm $SRC/bedtools-2.30.0.tar.gz && rm -r $SRC/bedtools2

## CTAT-VirusIntegrationFinder
WORKDIR $SRC
ENV VIF_CO="3e48b729c099827983c63b772823e5e4fc45f9cd"
RUN git clone --recursive https://github.com/broadinstitute/CTAT-VirusIntegrationFinder.git && \
        cd CTAT-VirusIntegrationFinder && \
        git checkout ${VIF_CO} && \
        git submodule init && git submodule update && \
        mv * $BIN && \
	    cd ${BIN}/util/bamsifter && ./build_htslib.sh && make

RUN wget -q https://github.com/broadinstitute/cromwell/releases/download/59/cromwell-59.jar -O $SRC/CTAT-VirusIntegrationFinder/WDL/cromwell-59.jar
RUN wget -q https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh -O $SRC/CTAT-VirusIntegrationFinder/WDL/monitor_script.sh
RUN chmod a+rx $SRC/CTAT-VirusIntegrationFinder/WDL/monitor_script.sh

RUN apt-get -qq -y remove automake build-essential cmake curl git gcc g++ wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log
