FROM ubuntu:18.04
MAINTAINER bhaas@broadinstitute.org
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
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

#~~~~~~~~~~~
# Python 
#~~~~~~~~~~~
RUN pip3 install --upgrade pip && pip3 install \
	igv-reports==1.0.2 \
	pandas==1.1.5 \
	argparse==1.4.0 

ENV SRC /usr/local/src
ENV BIN /usr/local/bin

#~~~~~~~~~~~
# R packages 
#~~~~~~~~~~~
RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("tidyverse")'
RUN R -e 'BiocManager::install("argparse")'
RUN R -e "install.packages(c('argparse', 'tidyr', 'tidyverse'), repos = 'http://cran.us.r-project.org')"

#~~~~~~~~~~~
## Samtools
#~~~~~~~~~~~
WORKDIR $SRC
ENV SAMTOOLS_VERSION=1.12
RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    cd $SRC && \
    wget $SAMTOOLS_URL && \
    tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    rm $SRC/samtools-${SAMTOOLS_VERSION}.tar.bz2 && rm -r $SRC/samtools-${SAMTOOLS_VERSION}

#~~~~~~~~~~~
## Picard tools
#~~~~~~~~~~~
WORKDIR $SRC
RUN wget https://github.com/broadinstitute/picard/releases/download/2.25.2/picard.jar
ENV PICARD_HOME $SRC/picard.jar

#~~~~~~~~~~~
## STAR Aligner
#~~~~~~~~~~~
WORKDIR $SRC
ENV STAR_VERSION=2.7.8a
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
    tar xf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
    mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin && \
    rm $SRC/${STAR_VERSION}.tar.gz && rm -r $SRC/STAR-${STAR_VERSION}

#~~~~~~~~~~~
## Bedtools
#~~~~~~~~~~~
WORKDIR $SRC
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
   tar xf bedtools-2.30.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/ && \
   rm $SRC/bedtools-2.30.0.tar.gz && rm -r $SRC/bedtools2


# blast
WORKDIR $SRC
ENV BLASTPLUS_VERSION 2.12.0
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLASTPLUS_VERSION}/ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    tar xvf ncbi-blast-${BLASTPLUS_VERSION}+-x64-linux.tar.gz && \
    cp ncbi-blast-${BLASTPLUS_VERSION}+/bin/* $BIN && \
    rm -r ncbi-blast-${BLASTPLUS_VERSION}+



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## CTAT-VirusIntegrationFinder
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WORKDIR $SRC
ENV CTAT_VIF_VER "1.3.0-predev"
ENV VIF_CO=e3439988e9f243b15ffad95325e82fcf889b7152
RUN git clone --recursive https://github.com/broadinstitute/CTAT-VirusIntegrationFinder.git && \
        cd CTAT-VirusIntegrationFinder && \
        git checkout ${VIF_CO} && \
        git submodule init && git submodule update && \
        make && \
        cp -r * $BIN 


#RUN wget -q https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh -O $SRC/CTAT-VirusIntegrationFinder/WDL/monitor_script.sh
#RUN chmod a+rx $SRC/CTAT-VirusIntegrationFinder/WDL/monitor_script.sh

RUN apt-get -qq -y remove automake build-essential cmake curl git gcc g++ wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

COPY sam_readname_cleaner.py $BIN/
