FROM bioconductor/bioconductor_docker

RUN echo "$GITHUB_PAT"

RUN --mount=type=secret,id=github_token \
  export GITHUB_PAT=$(cat /run/secrets/github_token) 
  
RUN echo "$GITHUB_PAT"

RUN apt-get update && apt-get -y upgrade && \
        apt-get install -y build-essential wget \
                libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev libcairo2-dev libxt-dev && \
        apt-get clean && apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN Rscript -e "install.packages('argparse')"
RUN Rscript -e "install.packages('R.utils')"
RUN Rscript -e "install.packages('magick')"
RUN Rscript -e "install.packages('devtools')"

RUN Rscript -e "library(devtools); install_github('shahcompbio/signals')"
RUN Rscript -e "library(devtools); install_github('caravagnalab/mobster')"

ADD policy.xml /etc/ImageMagick-6/policy.xml

#Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
        tar jxf samtools-1.9.tar.bz2 && \
        rm samtools-1.9.tar.bz2 && \
        cd samtools-1.9 && \
        ./configure --prefix $(pwd) && \
        make

ENV PATH=${PATH}:/usr/src/samtools-1.9

WORKDIR /usr/src