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


