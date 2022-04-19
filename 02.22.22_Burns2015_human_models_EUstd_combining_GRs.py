#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:40:56 2022

@author: adamo010
"""
#borrow heavily from 12.28.21_Muehlbauer2021_comm_GR_combining_files_and_adding_metadata.py

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

os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/EUstd_diet_model_outputs/")
    
filelist = glob.glob('*V01_Burns2015_human_EUstd_objective.txt')

sample_names_list = []
for file in filelist:
    sample_name = re.sub("_V01_Burns2015_human_EUstd_objective.txt", "", file) #delete extra crap on string.
    sep="_"
    sample_name2 = sample_name.split(sep, 1)[0] #need to remove extra crap on sample name string
    sample_names_list.append(sample_name2)
del file, sample_name, sample_name2, sep

df_list = []
for file in filelist:
    f = open(file, 'r')
    content_list2 = f.readlines()
    gr_num = str(content_list2[2])
    gr_num_list = gr_num.split(sep="= ")
    gr_num2 = round(float(gr_num_list[1]),3)   
    df_list.append(gr_num2)
del f, content_list2, gr_num, gr_num_list, gr_num2 

GR_dict = dict(zip(sample_names_list, df_list))

#nice. Save our output file
with open('02.22.22_Burns2015_human_models_EUstd_GRs_combined.csv', 'w') as f:
    for key in GR_dict.keys():
        f.write("%s,%s\n"%(key,GR_dict[key]))
del f, key, GR_dict       
        

