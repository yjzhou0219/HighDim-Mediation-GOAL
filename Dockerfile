FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    make \
    cmake \
    git \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /project

COPY R /project/R
WORKDIR /project/R

RUN Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(project=getwd(), prompt=FALSE)'
RUN Rscript -e 'library(lqa); library(DoubleML); library(ggh4x); packageVersion("lqa")'

CMD ["R"]
