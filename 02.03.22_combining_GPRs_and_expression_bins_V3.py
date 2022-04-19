#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 15:58:33 2022

@author: adamo010
"""

#01.20.22 update: trying to clean up my mess to make something sensible
#02.03.22 update: same.


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

##########step 1: import data
orig_GPR_file = pd.read_csv("RECON3_GPR_table_for_CORDA_full_V3.csv") #this is consistent. V3 is for typos.
#new addition: we have one Rx that is too large to have in a csv file, so it's in a text file. I'm going to try to add it manually.
text_file = open("Recon3_edited_ATPS4mi_Rx.txt", "r") #open text file in read mode
ATPS4mi_Rx = text_file.read() #read whole file to a string
text_file.close() #close file
del text_file
orig_GPR_file= orig_GPR_file.replace("HELPHELPHELP", ATPS4mi_Rx) #substitute the string ATPS4mi_Rx for "HELPHELPHELP"
#import all TPM abundance files
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto/")
RNAseq_file_list = []
RNAseq_file_names = []
for file in os.listdir():
    if str("expression_binned.csv") in str(file):
        RNAseq_file = pd.read_csv(file)
        RNAseq_file.drop(columns={"Unnamed: 0", "length", "eff_length", "est_counts", "tpm"}, inplace=True)
        RNAseq_file_list.append(RNAseq_file)
        RNAseq_file_names.append(str(file))
del RNAseq_file, file
RNAseq_file_names2 = [x.replace('_expression_binned.csv', '') for x in RNAseq_file_names]
RNAseq_dict = dict(zip(RNAseq_file_names2, RNAseq_file_list))
del RNAseq_file_names
#import table for matching geneids
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/")
genename_lookup_table_full = pd.read_csv("GRCh38_lookup_table.txt", delimiter = "\t")    
genename_lookup_table_abridged = genename_lookup_table_full[["name","geneName"]] #keep only relevant columns
del genename_lookup_table_full
#import table for ccds_ids
ccds_table_full = pd.read_csv("CCDS_lookup_table.txt", delimiter = "\t")
ccds_table_abridged = ccds_table_full[["gene", 'gene_id', 'ccds_id']]
ccds_table_abridged["ccds_id"] = ccds_table_abridged["ccds_id"].str.replace('CCDS','') #the CCDS ids have "CCDS" in front of them; that's not helpful.
del ccds_table_full

##############STEP 2 converting RNAseq data to something useable- basically, exchange weird gene ids into ccds
RNAseq_dict_ccds_str={} 

for key, value in RNAseq_dict.items():
    #first, attach gene ids to target_idsusing genename_lookup_table_abridged
    RNAseq_file_geneids = pd.merge(left=value, right= genename_lookup_table_abridged, how="left", left_on="target_id", right_on="name")
    RNAseq_file_geneids.drop(columns={"name"}, inplace=True)
    #then, add ccds_ids from gene names using ccds_table-abridged
    RNAseq_file_ccds = pd.merge(left=RNAseq_file_geneids, right=ccds_table_abridged, how="left", left_on="geneName", right_on="gene")
    RNAseq_file_ccds= RNAseq_file_ccds.drop(columns={"gene", "gene_id"})
    #then do some duplicate removal stuff
    RNAseq_lookup_key = RNAseq_file_ccds.drop(columns={"target_id", "geneName"}) #new, cleaner df
    RNAseq_lookup_key.drop_duplicates(subset= ["exp_score","ccds_id"], inplace=True) #drop cases where exp_score-ccds_id pairs are identical
    RNAseq_lookup_key = RNAseq_lookup_key.dropna(subset=['ccds_id']) #drop cases where ccds_id is nan (all -1/0/1 anyway)
    RNAseq_lookup_key_nodup = RNAseq_lookup_key.sort_values('exp_score', ascending=False).drop_duplicates('ccds_id').sort_index()
    #in cases where ccds_id has multiple options for a bin score (multiple transcripts per gene, I think), keep largest bin
    RNAseq_lookup_key_dict = RNAseq_lookup_key_nodup.set_index('ccds_id')['exp_score'].to_dict() #create dictionary
    keys_values = RNAseq_lookup_key_dict.items() #create placeholder dict_object
    RNAseq_lookup_key_dict_str = {str(key): str(value) for key, value in keys_values} #a dictionary of ccds_ids (keys) and their bins (values)
    RNAseq_dict_ccds_str.update({key: RNAseq_lookup_key_dict_str}) #append ccds_id-bin dictionary as a value to RNAseq_dict, with sampleid as the key
    #yes, it is a dictionary within a dictionary. 
del key, value, RNAseq_file_geneids, RNAseq_file_ccds,RNAseq_lookup_key,RNAseq_lookup_key_nodup,RNAseq_lookup_key_dict,keys_values, RNAseq_lookup_key_dict_str
#nice and fast with just one sample.     
#note: RNAseq_dict_ccds_str should have n key-value pairs, where n is number of samples in RNAseq_dict. 

##########Step 3: write and run a function for swapping ccds_ids with their bins
#the input for this function is the key-value pair from RNAseq_dict_ccds_str
def substituting_gprs(RNAseq_lookup_key_dict_str, sample_name):
    subrulecol = []
    sample_GPR_file = copy.deepcopy(orig_GPR_file) #orig_GPR_file is shared across all samples, so don't want to screw with it
    for row in sample_GPR_file["Rx_rule"]:
        if "(" in row: #if ( in row, that means a) there is a gene associated withthis Rx, and there are more than one genes associated
            rowlist= re.split("[,*]", row) #start by comma and asterisk-separating items; ",|*" means "split by , or *
            rowlist2 = []
            for item in rowlist:
                try: 
                    float(item) #if the item can be made into a number (i.e.isn't a bracket)
                    if item in RNAseq_lookup_key_dict_str.keys():
                        item2 = [val for key, val in RNAseq_lookup_key_dict_str.items() if item==key]
                        item3= item2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
                            #print(item2) 
                        rowlist2.append(item3)
                    else:
                        rowlist2.append("0") #if the ccds_id is not in our RNAseq data, it has unknown confidence- assign a value of 0, per instructions in CORDA paper
                except ValueError:
                    rowlist2.append(item) #this is where the items can't be converted to float (i.e. max/min)
            rowstring2 = ""
            rowstring2= ",".join([str(item3) for item3 in rowlist2]) #reassemble
            rowstring3 = rowstring2.replace("(,", "(")
            rowstring3= rowstring3.replace(",)",")")
            subrulecol.append(rowstring3)
        else: #if there is no ( in row, the GPR must be in a single gene (NOTE that gene-less GPRs are in a separate file and are dealt with later)
            if row in RNAseq_lookup_key_dict_str.keys():
                 row2 = [val for key, val in RNAseq_lookup_key_dict_str.items() if row==key] #return the bin associated with that single gene
                 row3 = row2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
                 subrulecol.append(row3)
            else:
                subrulecol.append("0") ##if the ccds_id is not in our RNAseq data, it has unknown confidence- assign a value of 0, per instructions in CORDA paper
    del row, rowlist, rowlist2, item, item2, item3, rowstring2, rowstring3, row2, row3              
    sample_GPR_file["edited_Rx_rule"]= subrulecol
    sample_GPR_file.to_csv("_GPRs_binned_V3.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned_V3.csv"):
        	dst=str(sample_name)+file
        	os.rename(src,dst)   
    return
##now, run it
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/")
for key, value in RNAseq_dict_ccds_str.items():
    substituting_gprs(value, key)

###########Step 4: adding non-gene GPRs
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/")
binned_file_list = []
binned_file_names = []
for file in os.listdir():
    if("_GPRs_binned_V3.csv" in file): #CHANGE AS NEEDED
        binned_file = pd.read_csv(file)
        binned_file.drop(columns={"Unnamed: 0", "Rx_rule"}, inplace=True)
        binned_file_list.append(binned_file)
        binned_file_names.append(str(file))
del binned_file, file

binned_file_names2 = [x.replace('_GPRs_binned_V3.csv', '') for x in binned_file_names] #CHANGE AS NEEDED
binned_dict = dict(zip(binned_file_names2, binned_file_list))
del binned_file_names
#import a dataframe of reactions without associated genes/GPRs. NOTE HERE: in V3, all of these have been assigned high confidence. The rationale
#is that these Rxs are part of the Recon3 model but have no associated genes, which means they must be essential for model running.
nongene_Rxs = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/RECON3_solved_GPRs_for_Rxs_without_genes_V3.csv") #USING V3
#all RIGHT, let's make another function and run this.
binned_dict2= copy.deepcopy(binned_dict)
#also, import a csv file 
for key, value in binned_dict2.items():
    Rx_solved = []
    for row in value["edited_Rx_rule"]:
        try:
            solvedval = (eval(row)) #THE MAGIC TRICK! evaluates max/min statements. only works if max/min has been filled out correctly.
            Rx_solved.append(solvedval)
        except SyntaxError: #raised if max/min has a ccid still hanging out, or is missing a comma or something
            Rx_solved.append(row)
        except TypeError: #raised if there's already a single value (e/g 0) that's being read as a string
            Rx_solved.append(row)
    value["solved_Rx_rule"]= Rx_solved
    value2= pd.concat([value, nongene_Rxs], axis=0)
    value2.to_csv("_GPRs_binned_solved_V3.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned_solved_V3.csv"):
        	dst=str(key)+file
        	os.rename(src,dst)   
del key, value, Rx_solved, row, solvedval, value2, file, dst, src
