
# Read RDS files into list
out = lapply(snakemake@input[["rds"]], function(RDS){
    # read in file
    readRDS(RDS)
})

# Set names as EFO IDs
names(out) = unlist(lapply(snakemake@input[["rds"]], function(RDS){
    stringr::str_remove(basename(RDS), ".rds")
}))

# Save
saveRDS(out, snakemake@output[[1]])