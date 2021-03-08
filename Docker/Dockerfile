FROM ubuntu:18.04

MAINTAINER bhaas@broadinstitute.org

RUN apt-get -qq update && apt-get -qq -y install automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
    fort77 \
    g++ \
    gcc \
    gfortran \
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
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-distutils \
    python3-setuptools \
    unzip \
    wget \
    x11-common \
    zlib1g \
    zlib1g-dev \
    zlibc


# make python3 be the default python
RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
  python get-pip.py

RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm install URI::Escape


ENV SRC /usr/local/src
ENV BIN /usr/local/bin
WORKDIR $SRC

## Samtools
ENV SAMTOOLS_VERSION=1.9
RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    cd $SRC && \
    wget $SAMTOOLS_URL && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install

## Bedtools
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
   tar -zxvf bedtools-2.30.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/ && \
   rm $SRC/bedtools-2.30.0.tar.gz

## STAR Aligner
ENV STAR_VERSION=2.7.2b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin

# IGV
ENV IGV_REPORTS_CO=1086ed7258c09bfe213836e2dc53b043c352a4a9
RUN git clone https://github.com/igvteam/igv-reports.git && \
    cd igv-reports && \
    git checkout ${IGV_REPORTS_CO} && \
    pip install -e .


# R
ENV DEBIAN_FRONTEND=noninteractive
ENV R_VERSION=R-3.6.3

RUN curl https://cran.r-project.org/src/base/R-3/$R_VERSION.tar.gz -o $R_VERSION.tar.gz && \
    tar xvf $R_VERSION.tar.gz && \
    cd $R_VERSION && \
    ./configure && make && make install

RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("tidyverse")'
RUN R -e 'BiocManager::install("argparse")'
RUN R -e ' BiocManager::install("tidyverse")'

COPY ./ /usr/local/bin/
RUN  cd /usr/local/bin/util/bamsifter && make

## Picard tools
RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.3/picard.jar
ENV PICARD_HOME $SRC

RUN pip install pandas

RUN wget https://github.com/billzt/gff3sort/archive/v1.0.0.zip && \
    unzip v1.0.0.zip && rm v1.0.0.zip && mv gff3sort-1.0.0/gff3sort.pl $BIN && mv gff3sort-1.0.0/Sort $BIN

RUN apt-get -qq -y remove build-essential curl gcc g++ git make unzip wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log
