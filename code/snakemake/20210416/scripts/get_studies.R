#!/usr/bin/env Rscript

# Libraries

library(gwasrapidd)

# Read in data

assocs = readRDS(snakemake@input[[1]])

# Get linked studies

studies_key = gwasrapidd::association_to_study(unique(assocs@associations$association_id))

# Get study data

studies = gwasrapidd::get_studies(study_id = unique(studies_key$study_id))

# Save to file

saveRDS(studies_key, snakemake@output[["key"]])
saveRDS(studies, snakemake@output[["studies"]])