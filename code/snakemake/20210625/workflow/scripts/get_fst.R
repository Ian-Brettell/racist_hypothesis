#!/usr/bin/env Rscript

# Libraries

library(tidyverse)
library(pegas)

# Read in genotypes

genos_raw = lapply(snakemake@input[["vcfs"]], function(VCF){
  n_variants = nrow(pegas::VCFloci(VCF))
  pegas::read.vcf(VCF, to = n_variants)
})
names(genos_raw) = basename(snakemake@input[["vcfs"]]) %>%
    stringr::str_remove(".vcf.gz")

# Read in population file
pop_file = readr::read_csv(snakemake@input[["pop_file"]])

# Create vector of populations
populations = unlist(lapply(rownames(genos_raw[[1]]), function(sample){
    pop_file$Population[pop_file$Sample == sample]
}))

# Calculate Fst

fst_pegas = lapply(genos_raw, function(TRAIT_ID){
  pegas::Fst(TRAIT_ID, pop = populations)
})

# Write to RDS file

dir.create(dirname(snakemake@output[[1]]), showWarnings = F, recursive = T)
saveRDS(fst_pegas, snakemake@output[[1]])