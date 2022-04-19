#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 10:27:10 2022

@author: adamo010
"""
#borrow heavily from 12.28.21_Muehlbauer2021_fluxes_combining_files_and_adding_metadata.py

import numpy as np
import cobra
import csv
import subprocess
import os
import pandas as pd
import fnmatch
import shutil
import glob
import re
import copy

######step 1: import files
metadata= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/CRC_metadata_rnaseq_edited_by_BA.csv")

os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/highfiber_diet_model_outputs/")
    
filelist = glob.glob('*V01_Burns2015_human_highfiber_fluxes.csv')

sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V01_Burns2015_human_highfiber_fluxes.csv", "", file) #delete extra crap on string.
    sep="_"
    sample_name2 = sample_name.split(sep, 1)[0] #need to remove extra crap on sample name string
    sample_names_list.append(sample_name2)
del file, sample_name, sample_name2, sep

df_list = []
for file in filelist:
    df= pd.read_csv(file)
    df_list.append(df)
del file, df

fluxpath_dict = dict(zip(sample_names_list, df_list))

#import metabolites key. NEW KEY. For peoples only.
#bring in metab key
fluxpath_key = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/recon_model_fluxpaths.tsv", sep="\t")
fluxpath_key_slimmed = fluxpath_key.loc[:, 'abbreviation':"subsystem"]
fluxpath_key_slimmed.drop(columns={'formula'}, inplace=True) #get a copy error here, but ignore it.
del fluxpath_key

#TESTING GROUNDS
test_file = df_list[0]
test_file = test_file.rename(columns={"Unnamed: 0": "fluxpath_name"})
#test_file= test_file.set_index('fluxpath_name') #set row names to flux pathway names  
test_file["sample_id"] = "s38"
test_file["fluxpath_name"] = test_file["fluxpath_name"].map(lambda x: x.replace("[c]", "").replace("[l]", "").replace("[r]", "").replace("[e]", ""))
test_file_plus_fluxpath_info = pd.merge(left=test_file, right=fluxpath_key_slimmed, how="left", left_on="fluxpath_name", right_on = "abbreviation")
test_file_plus_fluxpath_info.drop(columns={'abbreviation'}, inplace=True)
test_file_plus_metadata_too = pd.merge(left=test_file_plus_fluxpath_info, right=metadata, how="left",\
                                           left_on="sample_id", right_on = "Tissue_Tube_ID")
#nice.

modded_df_list = []
for key, value in fluxpath_dict .items():
    test_file = value
    test_file = test_file.rename(columns={"Unnamed: 0": "fluxpath_name"})
    test_file["sample_id"] = str(key)
    test_file["fluxpath_name"] = test_file["fluxpath_name"].map(lambda x: x.replace("[c]", "").replace("[l]", "").replace("[r]", "").replace("[e]", ""))
    test_file_plus_fluxpath_info = pd.merge(left=test_file, right=fluxpath_key_slimmed, how="left", left_on="fluxpath_name", right_on = "abbreviation")
    test_file_plus_fluxpath_info.drop(columns={'abbreviation'}, inplace=True)
    test_file_plus_metadata_too = pd.merge(left=test_file_plus_fluxpath_info, right=metadata, how="left",\
                                           left_on="sample_id", right_on = "Tissue_Tube_ID")
    modded_df_list.append(test_file_plus_metadata_too)
del key, value, test_file, test_file_plus_fluxpath_info, test_file_plus_metadata_too    

all_fluxpath_df = pd.concat(modded_df_list, axis = 0) #fuck it, we'll stack 'em and make it a longform df

all_fluxpath_df.to_csv("02.14.22_Burns2015_human_models_highfiber_fluxes_combined.csv")

