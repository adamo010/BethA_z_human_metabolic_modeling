#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 09:45:23 2022

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

#goal is to combine output from 02.09.22_setting_diets_for_human_models_allsamples.py; specifically, the growth rates
#starting point: 12.28.21_Muehlbauer2021_comm_GR_combining_files_and_adding_metadata.py

#first, import metadata
metadata= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/CRC_metadata_rnaseq_edited_by_BA.csv")

os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/EUstd_diet_model_outputs/")
    
GR_file_list = glob.glob('*V01_Burns2015_human_EUstd_objective.txt')

GR_file_names = []

for elem in GR_file_list:
    sample_name = re.sub("_V01_Burns2015_human_EUstd_objective.txt", "", elem)  
    GR_file_names.append(sample_name)        
del sample_name    

#testing one sample
test_GR_file = GR_file_list[0]
test_GR_file2 = pd.read_csv(test_GR_file) #produces a DF with two rows, one column. only want the first row
test_GR_string = test_GR_file2.iloc[0]['Objective'] #print the "single cell" to a string (yes, I know it's a bit of a farce)
sep = ' = ' #this is our delimiter: we want all text prior to this removed
test_GR_string2 = test_GR_string.split(sep, 1)[1] #split, and only keep text after this delimiter
test_GR_string3= round(float(test_GR_string2), 2) #convert to float and round to 2 decimal places
GRs.update({sample_name: test_GR_string3}) #this won't work now,but will in the forloop
#all right, single sample worked out. Now we make a dictionary.     
    
#START HERE
GRs =  {}
for elem in GR_file_list:
    sample_name = re.sub("_V01_Burns2015_human_EUstd_objective.txt", "", elem) #delete extra crap on string.
    GR_file = pd.read_csv(elem)
    GR_string = GR_file.iloc[1]['Objective'] #print the "single cell" to a string (yes, I know it's a bit of a farce)
    sep = ' = ' #this is our delimiter: we want all text prior to this removed
    GR_string2 = GR_string.split(sep, 1)[1] #split, and only keep text after this delimiter
    GR_string3= round(float(GR_string2), 2) #convert to float and round to 2 decimal places
    GRs.update({sample_name: GR_string3})
del sample_name, GR_file, GR_string, sep, GR_string2, GR_string3
    
#nice. Save our output file
with open('Burns2015_human_models_EUstd_GRs_combined.csv', 'w') as f:
    for key in GRs.keys():
        f.write("%s,%s\n"%(key,GRs[key]))
del f, key, GR_file_list, GRs, GR_file_names        
        
#excellent!