FROM rocker/tidyverse:4.0.0

# Install packages required for `sf`, which in turn is required for `pegas`
RUN apt update
RUN apt install -y libgdal-dev libudunits2-dev

# Install R packages
RUN install2.r --error \
    here \
    pegas \
    knitr \
    karyoploteR \
    ggridges \
    gwasrapidd \
    cowplot
