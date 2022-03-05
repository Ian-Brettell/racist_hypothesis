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

# Filtered traits (2045/3459 which had >= 1 SNP after clumping)
filt_traits = pd.read_table(config["filt_trait_ids_file"])
EFO_IDS_FILT = filt_traits["EFO_ID"]

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# Create combinations of clumping parameters and efo_ids
# Put into dictionary
in_dict = {'clump_r2': config["clump_r2"], 'clump_kb': config["clump_kb"]}
# Create data frame from dictionary
new_df = pd.DataFrame(in_dict)
# Repeat by length of EFO_IDS_FILT list
new_df_rep = pd.concat([new_df]*len(EFO_IDS_FILT), ignore_index = True)
# Add efo_ids column
new_df_rep["efo_id"] = EFO_IDS_FILT.values.tolist() * len(new_df)

# Create columns for zipping
EFO_IDS_FILT_ZIP = new_df_rep["efo_id"]
CLUMP_R2_ZIP = new_df_rep["clump_r2"]
CLUMP_KB_ZIP = new_df_rep["clump_kb"]
DATE_OF_COLLECTION_ZIP = [config["date_of_gwas_collection"]]*len(EFO_IDS_FILT_ZIP)




