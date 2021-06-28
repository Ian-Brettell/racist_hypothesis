#!/usr/bin/env Rscript

# Libraries

library(gwasrapidd)
library(tidyverse)

in_file = readRDS(snakemake@input[[1]])

out = dplyr::left_join(in_file@risk_alleles %>%
                            dplyr::select(association_id, variant_id),
                       in_file@associations %>%
                            dplyr::select(association_id, pvalue),
                       by = "association_id") %>%
        dplyr::select(SNP = variant_id, P = pvalue) %>%
        readr::write_tsv(snakemake@output[[1]])