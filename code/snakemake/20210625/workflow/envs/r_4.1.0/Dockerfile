FROM rocker/tidyverse:4.1.0

# Install packages required for R packages
RUN apt update
RUN apt install -y \
    # required for `sf`, which in turn is required for `pegas`
    libgdal-dev libudunits2-dev \
    # required for `pegas`
    libglpk40 \
    # to stop error in RStudio `libXt.so.6: cannot open shared object file: No such file or directory`
    libxtst6
    

# Install CRAN packages
RUN install2.r --error \
    here \
    pegas \
    knitr \
    ggridges \
    cowplot \
    remotes \
    plotly \
    DT

# Install `gwasrapidd`
RUN R -e "remotes::install_github('ramiromagno/gwasrapidd')"
