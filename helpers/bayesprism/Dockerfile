FROM rocker/tidyverse:4.3.0

# install libglpk-dev
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

# RUN install2.r BiocManager

RUN install2.r \
    BiocManager \
    arrow \
    fs \
    logger \
    optparse \
    readr

RUN Rscript -e "BiocManager::install('bluster')" \
    && Rscript -e "BiocManager::install('scran')" \
    && Rscript -e 'library("devtools"); install_github("Danko-Lab/BayesPrism/BayesPrism")'

COPY run_bayesprism.R /src/

ENTRYPOINT ["Rscript", "/src/run_bayesprism.R"]
