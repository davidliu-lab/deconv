FROM rocker/tidyverse:4.3.0

# install libglpk-dev
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libglpk-dev \
    && rm -rf /var/lib/apt/lists/*


RUN install2.r BiocManager

RUN Rscript -e "BiocManager::install('bluster')"

RUN Rscript -e "BiocManager::install('scran')"

RUN Rscript -e 'library("devtools"); install_github("Danko-Lab/BayesPrism/BayesPrism")'
#     && rm -rf /tmp/downloaded_packages/

RUN Rscript -e 'library("BayesPrism")'

RUN install2.r \\
    arrow \\
    fs \\
    logger \\
    optparse \\
    readr

COPY run_bayesprism.R /src/

ENTRYPOINT ["Rscript", "/src/run_bayesprism.R"]

# ENTRYPOINT ["Rscript", "/src/run_bayesprism.R"]
