#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 12:37:24 2021

@author: adamo010
"""
#so, in order to use CORDA to make human metabolic models, the expression of each gene first needs to be binned into
#high/medium/low/nonexpressed categories (there is also a zero category for genes not present in the dataset/model,
#but we'll deal with that when it inevitably arises). The goal here is to import the read count table that Sambhawa
#sent, see if there are natural breaks, and bin each gene appropriately. For now, we can do it for all samples,
#but we will likely have to split it up by sample to make sample-specific models. 

#also, fuck R

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

#first, switch to our NEW FOLDER
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto")

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

#add 1 to each value of TPM
for df in list_of_dataframes:
    #print(list(df.columns.values)) 
    df['tpm'] = df['tpm'] + 1

for df in list_of_dataframes:
    df['log_tpm'] = np.log10(df['tpm'])
    
for df in list_of_dataframes:
    df.to_csv("_adjusted.csv")
    
TPM_files_adj = dict(zip(list_of_samples, list_of_dataframes)) 

for key, value in TPM_files_adj.items():                      
    value.to_csv("_adjusted.csv")
    for file in os.listdir():                       
        src=file
        if fnmatch.fnmatch(file, "_adjusted.csv"):
            dst = str(key)+file
            os.rename(src,dst)

#graphing in python sucks. I'm exporting this and importing into R. 
    
#expression_data = pd.read_csv('all_samples_protein_coding_subread_counts.txt', delimiter = "\t")
#expression_data.rename(columns={"Unnamed: 0": "gene_id"}, inplace=True)
#
#now, plot data: can I make 10 plots of randomly selected sets of 100 rows?
#expression_data.plot.bar()
##########test sample: B02_S2
B02_S2= pd.read_csv("B02_S2_abundance.tsv", sep='\t')
B02_S2["tpm"] = B02_S2["tpm"] + 1
B02_S2['log_tpm'] = np.log10(B02_S2['tpm'])

#let's get graphing
graph = B02_S2.hist(column='log_tpm', bins = 250)

################### 03.24.21  updates #####################
#graphing in spyder sucks. But even when I graphed in R, nothing useful came out of it.
#After a stupid literature search, I have come up with the following bins:
#- Score of 0: genes not included in dataset. 
#- Score of -1: no expression: TPM>= 0.5 (this gels with a conservative take from A model based criterion for gene expression calls using RNA-seq data)
#- Score of 1: low expression: TPM between 0.5 and 10
#- Score of 2: medium expression: TPM between 11 and 1000 (wow, what a range)
#- Score of 3: high expression: over 1000 TPM    
    
#so, for ever single gene in every single sample, convert TPM values to 0/-1/1/2/3
#Note that this is done on untransformed data, so I'm re-importing the files

#first, switch to our NEW FOLDER
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto")

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


#oh, hang on, I can just add another column. exp_score.
def tpm_scoring(row):
    if row["tpm"] == 0:
        return 0
    elif 0 < row["tpm"] <= 0.5:
        return -1
    elif 0.5 < row["tpm"] <= 10:
        return 1
    elif 10 < row["tpm"] <= 1000:
        return 2
    elif row["tpm"] > 1000:
        return 3

#dataframe.apply((lambda row: tpm_scoring(row), axis=1))
dataframe['exp_score'] = dataframe.apply (lambda row: tpm_scoring(row), axis=1)
#neat. Works on one sample.

#let's run it for all of them.
for key, value in TPM_files.items():
    value['exp_score'] = value.apply (lambda row: tpm_scoring(row), axis=1)
   # value.to_csv("_expression_binned.csv")
    #for file in os.listdir():                       
        #src=file
       # if fnmatch.fnmatch(file, "_expression_binned.csv"):
           # dst = str(key)+file
            #os.rename(src,dst)

#now, get counts for everything
#test with one sample
counts_dict = dataframe['exp_score'].value_counts().to_dict()

list_of_count_dicts = []
list_of_dict_names = []

for key, value in TPM_files.items():
    count_dict = value['exp_score'].value_counts().to_dict()
    list_of_count_dicts.append(count_dict)
    list_of_dict_names.append(key)
    
dict_of_count_dicts = dict(zip(list_of_dict_names, list_of_count_dicts)) 

#convert list of dictionaries to a single dataframe.
test_df = pd.DataFrame(dict_of_count_dicts)
test_df.to_csv("03.24.21_counts_of_binned_RNAseq_data.csv")






