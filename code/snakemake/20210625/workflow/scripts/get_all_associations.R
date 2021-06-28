# Libraries

library(readr)
library(gwasrapidd)

# Read in all traits file

all_traits = readr::read_tsv(snakemake@input[[1]])

counter = 0
lapply(all_traits@traits$efo_id, function(EFO_ID) {
  counter <<- counter + 1
  # set output file name
  out_path = file.path(snakemake@params[["output_dir"]],
                       paste(EFO_ID, ".rds", sep = ""))
  # if output file doesn't already exist, get associations and save
  if (!file.exists(out_path)){
    out = gwasrapidd::get_associations(efo_id = EFO_ID)
    # save
    saveRDS(out, out_path)
  } else {
    print(paste(counter, ". ", EFO_ID, ": File already exists.", sep = ""))
  }

})