# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Libraries

library(tidyverse)

# Get variables

## Debug
all_mafs = "/nfs/research/birney/users/ian/hmn_fst/mafs/1kg/20201028/all/all.csv"
ref_snps = list("/nfs/research/birney/users/ian/hmn_fst/gwasrapidd/20210809/high_cov/plink/clumped/0.1_1000/1001114.clumped",
                "/nfs/research/birney/users/ian/hmn_fst/gwasrapidd/20210809/high_cov/plink/clumped/0.1_1000/EFO_0000096.clumped")
file_with_seeds = "code/snakemake/20210625/config/20210809_filtered_traits.txt"
percent_interval = 5
out_files = list("/nfs/research/birney/users/ian/hmn_fst/gwasrapidd/20210809/controls/references/0.1_1000/1001114.csv",
                 "/nfs/research/birney/users/ian/hmn_fst/gwasrapidd/20210809/controls/references/0.1_1000/EFO_0000096.csv")

## True
all_mafs = snakemake@input[["all_mafs"]]
ref_snps = snakemake@input[["ref_snps"]]
file_with_seeds = snakemake@input[["file_with_seeds"]]
percent_interval = as.numeric(snakemake@params[["percent_interval"]])
#seed = snakemake@params[["seed"]]
out_files = snakemake@output[["csvs"]]


# Read in files

## All 1KG snps
all_mafs_df = readr::read_csv(all_mafs,
                              col_types = c("cicccddddd"))
## Filter out...
all_mafs_df = all_mafs_df %>%
    # missing IDs
    dplyr::filter(ID != ".") %>%
    # multiple IDs
    dplyr::filter(!(grepl(";", ID)))

## Trait reference SNPs
ref_snps_list = purrr::map(ref_snps, function(FILE){
    readr::read_table(FILE,
                      col_types = "--c---------"
                      )
})
names(ref_snps_list) = ref_snps %>% 
    unlist() %>% 
    basename() %>% 
    stringr::str_remove(".clumped")

## File with random seeds
seeds = readr::read_tsv(file_with_seeds,
                        col_types = c("ci"))

# assign names to `out_files` to ensure correct matches
names(out_files) = out_files %>% 
    unlist() %>% 
    basename() %>% 
    stringr::str_remove(".csv")

# Bin by 5% MAF intervals based on AF_EUR

## Set breaks
breaks = seq(0, 1, by = (percent_interval / 100))
## Add column with bins
all_mafs_df = all_mafs_df %>%
    dplyr::mutate(AF_EUR_BIN = cut(AF_EUR, breaks = breaks, labels = F, include.lowest = T))

# Split by interval

all_mafs_list = split(all_mafs_df, f = all_mafs_df$AF_EUR_BIN)

# Pull out MAFs for ref SNPs

ref_mafs = purrr::map(ref_snps_list, function(EFO_ID){
    # Get SNP IDs for each trait
    target_snps = EFO_ID %>% 
        dplyr::pull(SNP)
    # Get MAFs for those SNPs
    all_mafs_df %>%
        dplyr::filter(ID %in% target_snps)
})

# Loop over each ID in `ref_mafs` and pull out random SNP from same bin

counter = 0
final_list = purrr::map(ref_mafs, function(EFO_ID){
    # set counter
    counter <<- counter + 1
    # get EFO_ID
    TARGET_EFO = names(ref_mafs)[counter]
    # set seed
    seed = seeds %>% 
        dplyr::filter(EFO_ID == TARGET_EFO) %>% 
        dplyr::pull(SEED)
    set.seed(seed)
    # Loop over each SNP
    purrr::map_dfr(1:nrow(ref_mafs[[counter]]), function(i){
        target_bin = as.character(ref_mafs[[counter]][i, "AF_EUR_BIN"])
        random_snp = all_mafs_list[[target_bin]] %>%
            dplyr::slice_sample(n = 1) %>%
            dplyr::select(CONTROL_CHROM = CHROM,
                          CONTROL_POS = POS,
                          CONTROL_ID = ID,
                          CONTROL_AF_EUR = AF_EUR,
                          CONTROL_AF_EUR_BIN = AF_EUR_BIN)
        # Bind into DF
        original = ref_mafs[[counter]][i, c("CHROM", "POS", "ID", "AF_EUR", "AF_EUR_BIN")]
        out = cbind(original, random_snp)
    
        return(out)
    }) 
})

# Write to csv
counter = 0
lapply(final_list, function(EFO_ID){
    counter <<- counter + 1
    # target EFO
    TARGET_EFO = names(final_list)[counter]
    # get output file name
    OUT_PATH = out_files[[TARGET_EFO]]
    
    readr::write_csv(EFO_ID, OUT_PATH)
})
