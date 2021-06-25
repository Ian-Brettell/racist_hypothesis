FROM rocker/tidyverse:4.1.0

# Install packages required for `sf`, which in turn is required for `pegas`
RUN apt update
RUN apt install -y libgdal-dev libudunits2-dev

# Install CRAN packages
RUN install2.r --error \
    here \
    pegas \
    knitr \
#    BiocManager \
    ggridges \
    cowplot \
    remotes
    
# Install `karyoploteR` (bioconductor)
#RUN R -e "BiocManager::install('karyoploteR')"

# Install `gwasrapidd`
RUN R -e "remotes::install_github('ramiromagno/gwasrapidd')"
