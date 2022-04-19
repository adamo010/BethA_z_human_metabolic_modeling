#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 17:00:54 2022

@author: adamo010
"""
#now that we've successfully (or at least attemptedly) binned the Burns RNAseq data into the appropriate categories, it's time to give
#building models the ol' college try. Again. Ten months later. MICOM was more fun, sue me.

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
from cobra.io import load_matlab_model
import cobra
from corda import CORDA

os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/")

#all right, our goal here is to combine our sample-specific expression bin files with the GPR list from RECON3. This will allow us to
#paramterize recon3 in a sample-specific way.

#our sample-specific expression bins are here: /Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto
#each file is named as samplename_expression_binned.csv

#also, we'll need the RECON3 GPR list; it is at Recon3_GPR_dict_CORDA_updated.csv

orig_GPR_file = pd.read_csv("Recon3_GPR_dict_CORDA_updated.csv", header=None)
orig_GPR_file.rename(columns={0: "Rx_name", 1: "Rx_rule"}, inplace=True)
#might have to edit the Rx name column

#at some point we'll want to do this iteratively, but for now, let's just do it with one test file.
test_file = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto/B20_S20_expression_binned.csv")
test_file.drop(columns={"Unnamed: 0"}, inplace=True)

#Oh BOY, now we get to match Rx name with target_id!!!! SO EXCITING!!!
expression_data = pd.read_csv('/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_analyzed_data/all_samples_protein_coding_subread_counts.txt', delimiter = "\t")
#let's try this table:
lookup_table_full = pd.read_csv("GRCh38_lookup_table.txt", delimiter = "\t")    
#all right, seems good. Need two columns here: name (ensemblid) and geneName
list(lookup_table_full.columns.values) #print column names
lookup_table_relevent = lookup_table_full[["name", "name2", "geneName", "geneName2"]]
#okay, probably just need name and geneName
lookup_table_relevent2 = lookup_table_full[["name","geneName"]]

#let's pair the lookup table with the testfile
test_file_geneids= pd.merge(left=test_file, right= lookup_table_relevent2, how="left", left_on="target_id", right_on="name")

#good show. now drop extra columns
test_file_geneids_clean = test_file_geneids.drop(columns={"length", "eff_length", "est_counts", "tpm", "name"})

#a new attempt to bring in a useful lookup table
ccds_table_full = pd.read_csv("CCDS_lookup_table.txt", delimiter = "\t")
ccds_table_abridged = ccds_table_full[["gene", 'gene_id', 'ccds_id']]
#the CCDS ids have "CCDS" in front of them; that's not helpful.
ccds_table_abridged["ccds_id"] = ccds_table_abridged["ccds_id"].str.replace('CCDS','')

#now, we need to merge ccds_table_abridged with test_file_geneids_clean
full_data_table = pd.merge(left=test_file_geneids_clean, right=ccds_table_abridged, how="left", left_on="geneName", right_on="gene")
full_data_table= full_data_table.drop(columns={"gene", "gene_id"})

#pfff, now for the fun part. We get to figure out how to solve equations in python.
#maybe... make a dictionary of ccds_id and exp_score? I don't think we can, because we'll end up with duplicates.
#can we remove duplicates first? the only thing we actually care about are the ccds_ids and exp_score, and those might be
#duplicated across unique target_ids.

test_lookup_key = full_data_table.drop(columns={"target_id", "geneName"})
test_lookup_key.drop_duplicates(subset= ["exp_score","ccds_id"], inplace=True)
#drop nans also- there are only three, that have exp_scores but no ccds_ids
test_lookup_key = test_lookup_key.dropna(subset=['ccds_id'])

#NEW: have to remove duplicate values from test_lookup_key
test_lookup_key2 = test_lookup_key.sort_values('exp_score', ascending=False).drop_duplicates('ccds_id').sort_index()

#convert to duct
test_lookup_key_dict = test_lookup_key2.set_index('ccds_id')['exp_score'].to_dict()

#we are now moving forward with test_lookup_key_dict and the orig_gpr file
#MMMMMMM, I feel like this might be easier to do in excel, huh?
#probably. BUT, I also don't want to have to do it for 88 samples. That'll take all day. 

#let's just see how long it will take to do ONE in excel.
test_lookup_key2.to_csv("B20_S20_scores_and_ccids.csv")
orig_GPR_file.to_csv("B20_S20_GPR_table.csv")
#the GPR table will be the same for all samples, but I want to save it with a sample-specific name so I can re-import with a sample-specific name

#test_lookup1 = orig_GPR_file.replace({"Rx_rule": test_lookup_key_dict}) #this did not work
#for key, value in test_lookup_key_dict.items():
    #for i in range(len(orig_GPR_file.Rx_rule)):
        #if key == orig_GPR_file.Rx_rule[i]:
            #print(value)

#orig_GPR_file.Rx_rule = orig_GPR_file.Rx_rule.map(test_lookup_key_dict) didn't work


#orig_GPR_file["Rx_rule"]=orig_GPR_file["Rx_rule"].apply(str)
keys_values = test_lookup_key_dict.items()
test_lookup_key_dict = {str(key): str(value) for key, value in keys_values}

test_GPR_file = copy.deepcopy(orig_GPR_file)

#for row in test_GPR_file["Rx_rule"]:
    #print(type(row))
    #for key, value in test_lookup_key_dict.items():
        #row = row.replace(value, key)
        #print(row)
#del row, key, value   
#this breaks a local run. So, let's make sure it works on a single row. 

#this might help also- delete all rows in test_GPR_file with nan values in Rx_rule
#ugh, there's only 48 actual good Rxs- this sucks
test_GPR_file = test_GPR_file[test_GPR_file["Rx_rule"].str.contains("nan")==False ]
 
testo = test_GPR_file[['Rx_rule']].iloc[0] #get one row as a tester string

#for row in test_GPR_file["Rx_rule"]:
    #print(type(row))
    #for key, value in test_lookup_key_dict.items():
        #row = row.replace(value, key)
        #print(row)
#del row, key, value   
#for key, value in test_lookup_key_dict.items():
    #testo = testo.replace(key, value)
    
#I think part of the problem is that I have a massive-ass dictionary to iterate over. I need to filter it by only ccds_ids present in
#my test_GPR_file   

for key, value in test_lookup_key_dict.items():
    testo= testo.replace({key: value}, regex=True)
#YESSSSSSSSS!!! THIS ONE WORKED!!!!!!!

#whew. okay. let's see how doing the whole df goes.

#test_GPR_file['Rx_rule_filled'] = test_GPR_file["Rx_rule"].replace(test_lookup_key_dict, regex=True)
#okay, that worked, but there seems to be a problem with the dictionary.... there are some very different values from 0/1/2/-3

#let's see if python can actually read those values as numbers
test_GPR_file["Rx_rule"]=test_GPR_file["Rx_rule"].apply(str)

#for row in test_GPR_file["Rx_rule"]:
    #res = [int(i) for i in row.split() if i.isdigit()]
    #print(res)

#hmmm. nothing is showing up as strings. I wonder if we need to use something other than replace...

#let's just... sanity check, make sure that the ccds ids in our test file are actually in the dictionary
#extract the column as a list, remove all instances of max/min/(), append each item to a BIG list, slim the list, and
#cross check with the dictionary

Rx_rule_list = test_GPR_file["Rx_rule"].tolist()
Rx_rule_list = [s.replace("max", "") for s in Rx_rule_list]
Rx_rule_list = [s.replace("min", "") for s in Rx_rule_list]
Rx_rule_list = [s.replace("(", "") for s in Rx_rule_list]
Rx_rule_list = [s.replace(")", "") for s in Rx_rule_list]

slimmed_Rx_rule_list = []
for item in Rx_rule_list:
    item2 = list(item.split(",")) #make a list out of the string, based on commas
    for elem in item2:
        if elem in slimmed_Rx_rule_list:
            pass
        else:
            slimmed_Rx_rule_list.append(elem)
del item, item2, elem

for item in slimmed_Rx_rule_list:
    if item in test_lookup_key_dict:
        pass
    else:
        print(item)

#OH GOOD THERE ARE LOTS OF CCDS IDS IN OUR TEST FILE THAT ARE NOT IN DICTIONARY. FFS.
#okay, okay. What does this actually mean? There are genes in our dataset that aren't in the model, right?
#because the test file is from RNAseq data and the dictionary is from the model.
#no, it's the opposite. slimmed_Rx_rule list is from the model, and the dictionary is from the sample.
#because there are some model items (slimmed_Rx_rule_list items) not present in the dictionary (sample), we're
#basically missing parameters.

#FUCK IT, let's call them zeroes. Will have to do this BEFORE any other subsequent iterative value subsetting,
#b3cause we only want one dictionary for the whole 88 samples. So we need missing values from all 88 samples
#to add to the dictionary.
#actually, no, we can probaly do it per sample, since the dictionary is per sample.

#create a list of absent ccds ids and match them to a (bin) value of 0
ccds_missing_from_test_sample = {}
for item in slimmed_Rx_rule_list:
    if item in test_lookup_key_dict:
        pass
    else:
        ccds_missing_from_test_sample[item]=0

#merge old dict from sample to new dict of 'missing' ccds
updated_test_lookup_key_dict = {**test_lookup_key_dict,**ccds_missing_from_test_sample}

#all right, let's try merging again.
test_GPR_file2 =  copy.deepcopy(orig_GPR_file)
test_GPR_file2 = test_GPR_file2[test_GPR_file2["Rx_rule"].str.contains("nan")==False ]
test_GPR_file2['Rx_rule_filled'] = test_GPR_file2["Rx_rule"].replace(updated_test_lookup_key_dict, regex=True)

#HHHHHHHHHHHHHHHHHHHHHHH okay that did not help. 
#let's try resub instead of replace
test_GPR_file3 =  copy.deepcopy(orig_GPR_file)
test_GPR_file3 = test_GPR_file3[test_GPR_file3["Rx_rule"].str.contains("nan")==False ]
test_GPR_file3['Rx_rule_filled'] = test_GPR_file3["Rx_rule"].replace(updated_test_lookup_key_dict, regex=True)

#convert all items in  updated_test_lookup_key_dict to strings first
keys_values = updated_test_lookup_key_dict.items()
updated_test_lookup_key_dict = {str(key): str(value) for key, value in keys_values}
del keys_values

for key, value in updated_test_lookup_key_dict.items():
    for row in test_GPR_file3["Rx_rule"]:
        expression = str(row)
        re.sub(key, value, expression)

#well, we've learned another way to get the same result. 

#let's see what these errors actually are. 
test_GPR_file3.to_csv("01.18.22_error_prone_expression_bin_filling.csv")

#so I think what's happening is that my lookup code is basically looking at "1576.1" and seeing "6.1", which becomes 2, which becomes 1576.1
#which means, really, it's a string issue. 

#I guess we'll have to expand out each row...
test_GPR_file4 =  copy.deepcopy(orig_GPR_file)
test_GPR_file4 = test_GPR_file4[test_GPR_file4["Rx_rule"].str.contains("nan")==False ]

bin_values = ["-1","0","1","2","3"]

for row in test_GPR_file4["Rx_rule"]:
    rowlist = row.split(sep=",") #start by comma-separating items
    #rowlist2 = []
    #for item in rowlist:
        #item2 = [val for key, val in updated_test_lookup_key_dict.items() if item in key]
        #rowlist2.append(item2)
    missing_items = []    
    rowlist2a = [updated_test_lookup_key_dict.get(item, item) for item in rowlist] #match values
    for item in rowlist2a:
        if item in bin_values:
            pass
        else:
            missing_items.append(item)     
    rowstring2a= ",".join(rowlist2a) #reassemble
    #rowlist3 = re.sep

#lame. So, rowlist2a/rowstring 2a works well. But it misses the first and last elements in each list. Need to match a substring
    
#additionally, there are several ccids that are PRESENT in the GPR file, and ABSENT in the RNAseq data. Meaning... this gene is not present in RNAseq
#therefore, these values should be 0 UNLESS they contain a max/min substring or a bracket. 

#something is wrong here. There are loads of ccds_ids that really, really don't exist. To the point that I'm not sure if they got edited weird.
#we might have to start a new file. This one is a mess. 












