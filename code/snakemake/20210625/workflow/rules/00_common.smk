#####################
# Libraries
#####################

import os.path
import pandas as pd

#####################
# Variables
#####################

# Config file
configfile: "code/snakemake/20210625/config/config.yaml"

CHRS  = config["contigs"]

# Date of GWAS data collection
DATE_OF_COLLECTION = config["date_of_gwas_collection"]

# Extensions
GZ_EXTS = config["gz_exts"]

## Only to run once, before anything else, to create file with traits and IDs
## Comment out the two substantive lines below, otherwise will throw an error
rule trait_ids_file:
    params:
        traits = config["traits"]
    output:
        config["trait_ids_file"]
    container:
        config["R"]
    script:
        "../scripts/trait_ids_file.R"

#rule get_all_traits:
#    output:
#        config["all_trait_ids_file"]
#    container:
#        config["R"]
#    script:
#        "../scripts/get_all_traits.R"

# Traits
traits = pd.read_table(config["trait_ids_file"])
EFO_IDS = traits['efo_id']

filtered_traits = pd.read_table(config["filt_trait_ids_file"])
EFO_IDS_FILT = filtered_traits['efo_id']

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()


