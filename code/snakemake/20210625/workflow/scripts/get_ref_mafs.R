#!/usr/bin/env Rscript

# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Collect variables

all_mafs = snakemake@input[["all_mafs"]]
ref_snps = snakemake@input[["ref_snps"]]
percent_interval = as.numeric(snakemake@params[["percent_interval"]])
out_file = snakemake@output[[1]]

# Libraries

library(tidyverse)

# Read in files

all_mafs_df = readr::read_csv(all_mafs,
                              col_types = c("cicccddddd"))

ref_snps_df = readr::read_delim(ref_snps,
                                delim = " ",
                                col_types = "--c---------",
                                trim_ws = T)

# Bin by 5% MAF intervals based on EUR_AF
#test = all_mafs_df %>% dplyr::slice_sample(n = 10000)

## Set breaks
breaks = seq(0, 1, by = (percent_interval / 100))
## Add column with bins
all_mafs_df = all_mafs_df %>%
    dplyr::mutate(EUR_AF_BIN = cut(EUR_AF, breaks = breaks, labels = F, include.lowest = T))

# Split by interval

all_mafs_list = split(all_mafs_df, f = all_mafs_df$EUR_AF_BIN)

# Pull out MAFs for ref SNPs

ref_mafs = all_mafs_df %>%
    dplyr::filter(ID %in% ref_snps_df$SNP)

# Loop over each ID in `ref_mafs` and pull out random SNP from same bin

set.seed(10)
purrr::map_dfr(1:nrow(ref_mafs), function(i){
    target_bin = as.character(ref_mafs[i, "EUR_AF_BIN"])
    random_snp = all_mafs_list[[target_bin]] %>%
        dplyr::slice_sample(n = 1) %>%
        dplyr::select(CONTROL_CHROM = CHROM,
                      CONTROL_POS = POS,
                      CONTROL_ID = ID,
                      CONTROL_EUR_AF = EUR_AF,
                      CONTROL_EUR_AF_BIN = EUR_AF_BIN)
    # Bind into DF
    original = ref_mafs[i, c("CHROM", "POS", "ID", "EUR_AF", "EUR_AF_BIN")]
    out = cbind(original, random_snp)

    return(out)
}) %>%
    readr::write_csv(out_file)

save.image()