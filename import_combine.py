################################################################################
#   Script that imports taxonomic data and metadata, merges the two and creates
#   subsets for UC, CD and control
#
#   inputs:
#       - meta_file: csv file containing meta data
#       - tax_file: tsv file containing the taxonomic profiles
#
#   outputs:
#       - OTUID_taxonomy_dic: reference from OTUID to bacteria names
#       - combined_data.csv: merged data file
#       - data_CD: contains CD data
#       - data_UC: contains UC data
#       - data_cont: contains data of controls
#
################################################################################

import numpy as np
import pandas as pd

meta_file = "../01 Raw data/hmp2_metadata_2018-08-20.csv"
tax_file = "../01 Raw data/taxonomic_profiles.tsv"

meta = pd.read_csv(meta_file, sep=',')
tax = pd.read_csv(tax_file, sep='\t')

# remove column with project columns
meta.drop(columns = "Project", inplace = True)

# tax.drop(columns = "#OTU ID", inplace = True)
OTUID_taxonomy_dic = tax[["#OTU ID","taxonomy"]]
OTUID_taxonomy_dic.to_pickle("01 temp data/OTUID_taxonomy_dic")

# remove taxonomy column
tax.drop(columns = "taxonomy", inplace = True)

# only want S16 data
meta_S16 = meta[meta["data_type"] == "biopsy_16S"]

# transpose taxonomic data so that IDs are in one column
tax = tax.transpose()


# use bacteria descriptions as the headers & drop last row with the descriptions
tax.columns = tax.iloc[0]
tax = tax[1:]

# insert one column --> assign indices (IDs)
tax.insert(0,"External ID","NAN")
tax["External ID"] = tax.index

comb_data = pd.merge(meta_S16,tax, on="External ID")

# drop all empty columns
comb_data.dropna(axis = 1, how = "all", inplace = True)

# save as csv
comb_data.to_csv("01 temp data/combined_data.csv")
# save as pickle
comb_data.to_pickle("01 temp data/combined_data")
