library(gwasrapidd)

# Get traits

all_traits = gwasrapidd::get_traits()

# Write to file

readr::write_tsv(all_traits@traits,
                 snakemake@output[[1]])