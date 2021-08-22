FROM bioconductor/bioconductor_docker

RUN Rscript -e "install.packages(pkgs = c('tidyverse', \
                                          'RColorBrewer', \
                                          'data.table', \
                                          'cowplot', \
                                          'rmarkdown', \
                                          'uwot', \
                                          'BiocManager'), \
                                        repos='https://cran.revolutionanalytics.com/', \
                                        dependencies=TRUE, \
                                        clean = TRUE)"

RUN apt-get update && apt-get -y upgrade && \
        apt-get install -y build-essential wget \
                libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev libcairo2-dev && \
        apt-get clean && apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN Rscript -e "install.packages('argparse')"
RUN Rscript -e "install.packages('R.utils')"
RUN Rscript -e "install.packages('magick')"

RUN Rscript -e "BiocManager::install('ggtree')"
RUN Rscript -e "BiocManager::install('ComplexHeatmap')"
RUN Rscript -e "BiocManager::install('IRanges')"
RUN Rscript -e "BiocManager::install('GenomicRanges')"

RUN Rscript -e "library(devtools); install_github('caravagn/pio')"
RUN Rscript -e "library(devtools); install_github('caravagn/easypar')"
RUN Rscript -e "library(devtools); install_github('caravagnalab/mobster')"
RUN Rscript -e "library(devtools); install_github('caravagnalab/VIBER')"
RUN Rscript -e "library(devtools); install_github('VPetukhov/ggrastr')"
RUN Rscript -e "library(devtools); install_github('shahcompbio/schnapps')"

WORKDIR /usr/src

#Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
        tar jxf samtools-1.9.tar.bz2 && \
        rm samtools-1.9.tar.bz2 && \
        cd samtools-1.9 && \
        ./configure --prefix $(pwd) && \
        make

ENV PATH=${PATH}:/usr/src/samtools-1.9