#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 18:57:10 2022

@author: adamo010
"""
import scipy
import numpy as np
import cobra
import csv
import subprocess
import os
import pandas as pd
import fnmatch
import matplotlib.pyplot as plt
import shutil
import glob
import re
import copy
import pickle

#goal is to combine output from 02.09.22_setting_diets_for_human_models_allsamples.py; specifically, the uptakes and secretions
#starting point: 12.28.21_Muehlbauer2021_indiv_spp_GRs_combining_files_and_adding_metadata.py

#first, import metadata
metadata= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/CRC_metadata_rnaseq_edited_by_BA.csv")
#EDIT ME TO MATCH sample_names_list !!!!!!!!!!!!!!!!!

os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/highfiber_diet_model_outputs/")
    
filelist = glob.glob('*V01_Burns2015_human_highfiber_secretion.csv')

sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V01_Burns2015_human_highfiber_secretion.csv", "", file) #delete extra crap on string.  
    sample_names_list.append(sample_name)
del file, sample_name

#NOTE for sample_names_list: need to remove all characters after the underscore.
sep = "_"
sample_names_list2 =[]
for item in sample_names_list:
    item2 = item.split(sep, 1)[0]
    sample_names_list2.append(item2)
del sep, sample_names_list, item, item2    

df_list = []
for file in filelist:
    df= pd.read_csv(file)
    df_list.append(df)
del file, df

secretion_dict = dict(zip(sample_names_list2, df_list))

#TEST HERE
#test_file = df_list[0]
#test_file = test_file.drop(columns=["Unnamed: 0"]) #delete extra column
#test_file["sample_id"] ="s40" #add new column with sampleid
#NEW: remove [c], [l], [r] from Metabolite column strings. 
#test_file["Metabolite"] = test_file["Metabolite"].map(lambda x: x.replace("[c]", "").replace("[l]", "").replace("[r]", "").replace("[e]", ""))

#bring in metab key
metab_key = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/recon_model_metabolites.tsv", sep="\t")
metab_key_slimmed = metab_key.loc[:, 'abbreviation':'fullName']    #use : to select full column
del metab_key

#test_file_plus_metab_info = pd.merge(left=test_file, right=metab_key_slimmed, how="left", left_on="Metabolite", right_on = "abbreviation")

#Now try to merge in metadata also
#test_file_plus_metadata_too = pd.merge(left=test_file_plus_metab_info, right=metadata, how="left", left_on="sample_id", right_on = "Tissue_Tube_ID")
#thumbs up. 

#REAL CODE
modded_df_list = []
for key, value in secretion_dict.items():
    test_file = value
    test_file = test_file.drop(columns=["Unnamed: 0"]) #delete extra column
    test_file["sample_id"] = str(key)
    test_file["Metabolite"] = test_file["Metabolite"].map(lambda x: x.replace("[c]", "").replace("[l]", "").replace("[r]", "").replace("[e]", ""))
    test_file_plus_metab_info = pd.merge(left=test_file, right=metab_key_slimmed, how="left", left_on="Metabolite", right_on = "abbreviation")
    test_file_plus_metadata_too = pd.merge(left=test_file_plus_metab_info, right=metadata, how="left",\
                                           left_on="sample_id", right_on = "Tissue_Tube_ID")
    modded_df_list.append(test_file_plus_metadata_too)
del key, value, test_file, test_file_plus_metab_info, test_file_plus_metadata_too    

all_secretion_df = pd.concat(modded_df_list, axis = 0) #fuck it, we'll stack 'em and make it a longform df

all_secretion_df.to_csv("02.14.22_Burns2015_human_models_highfiber_secretions_combined.csv")
    
##################and again for uptake!
filelist = glob.glob('*V01_Burns2015_human_highfiber_uptake.csv')

sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V01_Burns2015_human_highfiber_uptake.csv", "", file) #delete extra crap on string.  
    sample_names_list.append(sample_name)
del file, sample_name

sep = "_"
sample_names_list2 =[]
for item in sample_names_list:
    item2 = item.split(sep, 1)[0]
    sample_names_list2.append(item2)
del sep, sample_names_list, item, item2    

df_list = []
for file in filelist:
    df= pd.read_csv(file)
    df_list.append(df)
del file, df

uptake_dict = dict(zip(sample_names_list2, df_list))

#bring in metab key
metab_key = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/recon_model_metabolites.tsv", sep="\t")
metab_key_slimmed = metab_key.loc[:, 'abbreviation':'fullName']    #use : to select full column
del metab_key

#REAL CODE
modded_df_list = []
for key, value in uptake_dict.items():
    test_file = value
    test_file = test_file.drop(columns=["Unnamed: 0"]) #delete extra column
    test_file["sample_id"] = str(key)
    test_file["Metabolite"] = test_file["Metabolite"].map(lambda x: x.replace("[c]", "").replace("[l]", "").replace("[r]", "").replace("[e]", ""))
    test_file_plus_metab_info = pd.merge(left=test_file, right=metab_key_slimmed, how="left", left_on="Metabolite", right_on = "abbreviation")
    test_file_plus_metadata_too = pd.merge(left=test_file_plus_metab_info, right=metadata, how="left",\
                                           left_on="sample_id", right_on = "Tissue_Tube_ID")
    modded_df_list.append(test_file_plus_metadata_too)
del key, value, test_file, test_file_plus_metab_info, test_file_plus_metadata_too    

all_uptake_df = pd.concat(modded_df_list, axis = 0) #fuck it, we'll stack 'em and make it a longform df

all_uptake_df.to_csv("02.14.22_Burns2015_human_models_highfiber_uptakes_combined.csv")
    
    
    
    











