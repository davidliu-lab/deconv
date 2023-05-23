sudo apt install r-base r-base-dev

# create renv in current working directory
Rscript -e 'install.packages("renv", repos="http://cran.us.r-project.org")'
Rscript -e 'renv::init()'

# install packages
Rscript -e 'install.packages("littler", repos="http://cran.us.r-project.org")'
install2.r devtools
# Rscript -e 'install.packages("devtools", repos="http://cran.us.r-project.org")'
install2.r tidyverse
# Rscript -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")'
install2.r BiocManager
RUN Rscript -e "BiocManager::install('bluster')"
RUN Rscript -e "BiocManager::install('scran')"
RUN Rscript -e 'library("devtools"); install_github("Danko-Lab/BayesPrism/BayesPrism")'
