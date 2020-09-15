FROM ubuntu:18.04


MAINTAINER bhaas@broadinstitute.org


RUN apt-get update && apt-get install -y gcc g++ perl python3 automake make \
                                       wget curl libdb-dev \
                                       bzip2 zlibc zlib1g zlib1g-dev  default-jre \
                       python3-setuptools python3-dev build-essential python3-distutils \
                                       unzip libbz2-dev  liblzma-dev && \
    apt-get clean


# make python3 be the default python
RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
  python get-pip.py

RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm install URI::Escape

ENV SRC /usr/local/src
ENV BIN  /usr/local/bin

## Samtools
ENV SAMTOOLS_VERSION=1.7

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
   cd $SRC && \
      wget $SAMTOOLS_URL && \
         tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
            cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
               cd ../ && ./configure --without-curses && make && make install
               

## STAR Aligner
WORKDIR $SRC

ENV STAR_VERSION=2.7.2b
RUN STAR_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz" &&\
    wget -P $SRC $STAR_URL &&\
        tar -xvf $SRC/${STAR_VERSION}.tar.gz -C $SRC && \
            mv $SRC/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /usr/local/bin


## update igv-reports to current bleeding edge
RUN apt-get update && apt-get install -y git && apt-get clean

WORKDIR $SRC
ENV IGV_REPORTS_CO=1086ed7258c09bfe213836e2dc53b043c352a4a9
RUN git clone https://github.com/igvteam/igv-reports.git && \
    cd igv-reports && \
    git checkout ${IGV_REPORTS_CO} && \
    pip install -e .


#####
# Install R

WORKDIR $SRC

ENV DEBIAN_FRONTEND=noninteractive

ENV R_VERSION=R-3.6.3

RUN apt-get install -y gfortran fort77 libreadline-dev libcurl4-openssl-dev libx11-dev libxt-dev x11-common libcairo2-dev libpng-dev libreadline-dev libjpeg-dev pkg-config libtbb-dev  cmake rsync libssl-dev tzdata

RUN curl https://cran.r-project.org/src/base/R-3/$R_VERSION.tar.gz -o $R_VERSION.tar.gz && \
        tar xvf $R_VERSION.tar.gz && \
        cd $R_VERSION && \
        ./configure && make && make install


RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("tidyverse")'


## CTAT-VIF

WORKDIR $SRC

COPY ./CTAT-VirusIntegrationFinder/  /usr/local/bin/


## Samtools
ENV SAMTOOLS_VERSION=1.9

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
   cd $SRC && \
      wget $SAMTOOLS_URL && \
         tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
            cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
               cd ../ && ./configure --without-curses && make && make install


RUN R -e 'BiocManager::install("argparse")'


RUN apt-get install -y  libxml2-dev 
RUN R -e ' BiocManager::install("tidyverse")'
