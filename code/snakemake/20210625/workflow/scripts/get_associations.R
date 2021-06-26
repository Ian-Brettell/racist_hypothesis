#!/usr/bin/env Rscript

# Load libraries

library(gwasrapidd)

# Get traits

gcat_assoc = gwasrapidd::get_associations(efo_id = snakemake@params[["efo_id"]])

# Get list of associations

out = saveRDS(gcat_assoc, snakemake@output[[1]])
