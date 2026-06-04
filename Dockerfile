FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    gfortran make pandoc \
    git\
    && rm -rf /var/lib/apt/lists/*

WORKDIR /project

COPY R /project/R
WORKDIR /project/R

RUN Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(project=getwd(), prompt=FALSE)'
RUN Rscript -e 'library(lqa); packageVersion("lqa")'

CMD ["R"]