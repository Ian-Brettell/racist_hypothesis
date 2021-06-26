#!/usr/bin/env Rscript

library(gwasrapidd)
library(readr)
 
# Get traits and IDs
out = gwasrapidd::get_traits(efo_trait = snakemake@params[["traits"]]) 

# Extract data frame and write to file
out@traits %>% 
  readr::write_tsv(snakemake@output[[1]])