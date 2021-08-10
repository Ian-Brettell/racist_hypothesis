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

# Traits
traits = pd.read_table(config["all_trait_ids_file"])
EFO_IDS = traits['efo_id']

filtered_traits = pd.read_table(config["filt_trait_ids_file"])
EFO_IDS_FILT = filtered_traits['efo_id']

# Get additional traits (that had fewer than 50 unique SNP associations)
additional_traits = traits.loc[traits['trait'].isin(config["additional_traits"]), 'efo_id'].tolist()

# Append to list of filtered traits
EFO_IDS_FILT_PLUS = EFO_IDS_FILT.tolist() + additional_traits

# Add control
EFO_IDS_FILT_PLUS_CONTROL = EFO_IDS_FILT_PLUS + ["CONTROL"]

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()


