# lightweight docker image with latest version of R, several R dependencies
# and bayesprism package installed

FROM rocker/r-ver:4.3.0

RUN installGithub.r Danko-Lab/BayesPrism \\
    && rm -rf /tmp/downloaded_packages/

# entrypoint is a R session with bayesprism package loaded

ENTRYPOINT ["R", "-e", "library(bayesprism)"]

