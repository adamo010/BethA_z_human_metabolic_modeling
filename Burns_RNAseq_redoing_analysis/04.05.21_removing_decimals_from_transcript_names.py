#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 16:32:40 2021

@author: adamo010
"""
import os
import corda
import numpy as np
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
import scipy

#the goal is to remove all decimal places from the .tsv abundance files that were outputted from my kallisto rerun of the Burns data
#these decimals refer to the transcript version and aren't necessary for anything (as far as I can tell.)

#this might also be a nice opportunity to rename the files according to their sample names instead of.. the B25_S25 crap.

#set wd.
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto")

#import files
list_of_samples = []
list_of_dataframes = []

for file in os.listdir():
    if file.endswith(".tsv"): #there's an extra thing in this directory that I can't find or delete
        filename= str(file)
        samplename = re.sub("_abundance.tsv", "", filename)
        list_of_samples.append(samplename)
        dataframe = pd.read_csv(file, sep='\t')
        list_of_dataframes.append(dataframe)

TPM_files = dict(zip(list_of_samples, list_of_dataframes)) 
del file, filename, samplename    

#first, generate a list of sample names without the underscores.
list_of_samples_abv = []
for item in list_of_samples:
    split_string = item.split("_", 1)
    substring = split_string[0]
    list_of_samples_abv.append(substring)
del item, split_string, substring    

TPM_files2 = dict(zip(list_of_samples_abv, list_of_dataframes)) 

for key, value in TPM_files2.items():
    value["target_id"] = value["target_id"].str.replace(r'.', '')
    #print(value)
    value.to_csv("_transcript_trimmed.tsv", index=False, sep="\t")
    for file in os.listdir():                       
        src=file
        if fnmatch.fnmatch(file, "_transcript_trimmed.tsv"):
            dst = str(key)+file
            os.rename(src,dst)

#nice     
#04.06.21: had to add "index=False" to exported csvs, because it was confusing tximport to have an index column
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
