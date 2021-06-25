FROM rocker/tidyverse:4.1.0

# Install libgdal for `sf`, which in turn is needed for `pegas`
RUN apt install gdal-bin libgdal-dev

# Install R packages
RUN install2.r --error \
    here \
    pegas \
    knitr \
    karyoploteR \
    ggridges \
    gwasrapidd \
    cowplot
