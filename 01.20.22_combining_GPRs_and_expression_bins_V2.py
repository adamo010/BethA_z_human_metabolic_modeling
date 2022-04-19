#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 17:02:33 2022

@author: adamo010
"""
#01.20.22 update: trying to clean up my mess to make something sensible

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
#RECON3 GPR list; it is at Recon3_GPR_dict_CORDA_updated.csv
orig_GPR_file = pd.read_csv("RECON3_GPR_table_for_CORDA_full.csv") #NEW FILE!!
#orig_GPR_file.rename(columns={0: "Rx_name", 1: "Rx_rule"}, inplace=True) #no longer needed
#this is our BASE source file; each sample will get different parameters assigned to this file.
#the Rx_rule column contains several ccds_ids, which need to be exchanged for sample-specific -1/0/1/2/3 values. 
#if a ccds_id appears in the orig_GPR_file but is NOT in the RNA-seq data, the value should be 0
#if a ccds_id appears in the RNA-seq data but is NOT in the orig_GPR_file, drop that gene from RNA-seq data.

#at some point we'll want to do this iteratively, but for now, let's just do it with one test file.
test_RNAseq_file = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto/B20_S20_expression_binned.csv")
test_RNAseq_file.drop(columns={"Unnamed: 0", "length", "eff_length", "est_counts", "tpm"}, inplace=True)

###########step 2: converting RNAseq data to something useable
#this has to be done in two steps: match ensembleid to geneid, then geneid to ccds_id
#first, import table for matching geneids
genename_lookup_table_full = pd.read_csv("GRCh38_lookup_table.txt", delimiter = "\t")    
#keep only relevant columns
genename_lookup_table_abridged = genename_lookup_table_full[["name","geneName"]]
del genename_lookup_table_full
#merge with test_RNAseq_file
test_RNAseq_file_geneids= pd.merge(left=test_RNAseq_file, right= genename_lookup_table_abridged, how="left", left_on="target_id", right_on="name")
test_RNAseq_file_geneids.drop(columns={"name"}, inplace=True)
del test_RNAseq_file, genename_lookup_table_abridged

#second, import table for ccds_ids
ccds_table_full = pd.read_csv("CCDS_lookup_table.txt", delimiter = "\t")
ccds_table_abridged = ccds_table_full[["gene", 'gene_id', 'ccds_id']]
#the CCDS ids have "CCDS" in front of them; that's not helpful.
ccds_table_abridged["ccds_id"] = ccds_table_abridged["ccds_id"].str.replace('CCDS','')
del ccds_table_full

#step 2: merge ccds table with test_file (full_data_table --> test_RNAseq_file_ccds)
test_RNAseq_file_ccds = pd.merge(left=test_RNAseq_file_geneids, right=ccds_table_abridged, how="left", left_on="geneName", right_on="gene")
test_RNAseq_file_ccds= test_RNAseq_file_ccds.drop(columns={"gene", "gene_id"})
del test_RNAseq_file_geneids, ccds_table_abridged

##########step 3: create a lookup dictionary of ccds_ids and exp_score
test_RNAseq_lookup_key = test_RNAseq_file_ccds.drop(columns={"target_id", "geneName"})
test_RNAseq_lookup_key.drop_duplicates(subset= ["exp_score","ccds_id"], inplace=True)
#drop nans also- there are only three, that have exp_scores but no ccds_ids
test_RNAseq_lookup_key = test_RNAseq_lookup_key.dropna(subset=['ccds_id'])

#NEW: have to remove duplicate values from test_lookup_key
test_RNAseq_lookup_key_nodup = test_RNAseq_lookup_key.sort_values('exp_score', ascending=False).drop_duplicates('ccds_id').sort_index()
del test_RNAseq_lookup_key

#convert to dict
test_RNAseq_lookup_key_dict = test_RNAseq_lookup_key_nodup.set_index('ccds_id')['exp_score'].to_dict()
del test_RNAseq_lookup_key_nodup

#I'll also create a second dict that has only strings
keys_values = test_RNAseq_lookup_key_dict.items()
test_RNAseq_lookup_key_dict_str = {str(key): str(value) for key, value in keys_values}
del keys_values

###########step 4: testing.... how to exchange the ccds_ids in orig_GPR_file for their corresponding exp_scores in test_RNAseq_lookup_key_dict
#I guess we'll have to expand out each row...
test_GPR_file1 = copy.deepcopy(orig_GPR_file)
test_GPR_file1 = test_GPR_file1.head(10) #just take the first 10 rows, to speed up testing
#test_GPR_file1 = test_GPR_file1[test_GPR_file1["Rx_rule"].str.contains("nan")==False ] #delete nans- shouldnt' be any now.
#okay, I have a lot of questions about this file. we should have more than 49 Rxs with GPRs. 
#FIXED! with new file names

#bin_values = ["-1","0","1","2","3"]

subrulecol = []
for row in test_GPR_file1["Rx_rule"]:
    if "(" in row:
        #print(row)
        rowlist= re.split("[,*]", row) #start by comma and asterisk-separating items; ",|*" means "split by , or *
        #for item in rowlist:
            #item2 = re.split("*", item)
            #print(item2)
        #print(rowlist)
        rowlist2 = []
        for item in rowlist:
            try: 
                float(item)
                if item in test_RNAseq_lookup_key_dict_str.keys():
                    item2 = [val for key, val in test_RNAseq_lookup_key_dict_str.items() if item==key]
                    item3= item2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
                    #print(item2) 
                    rowlist2.append(item3)
                else:
                    rowlist2.append("0")
            except ValueError:
                    rowlist2.append(item) #this is where the items can't be converted to float
        rowstring2 = ""
        #print(rowlist2)
        rowstring2= ",".join([str(item3) for item3 in rowlist2]) #reassemble
        rowstring3 = rowstring2.replace("(,", "(")
        rowstring3= rowstring3.replace(",)",")")
        #print(rowstring2)
        #print(rowstring3)
        #print("     ")
        subrulecol.append(rowstring3)
    else:
        if row in test_RNAseq_lookup_key_dict_str.keys():
             row2 = [val for key, val in test_RNAseq_lookup_key_dict_str.items() if row==key]
             row3= row2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
             subrulecol.append(row3)
        else:
            subrulecol.append(row)
del row, rowlist, rowlist2, item, item2, item3, rowstring2, rowstring3, row2, row3              

test_GPR_file1["edited_Rx_rule"]= subrulecol

#cool. Now, let's make this into a function SEE BELOW...
######################################################################################
#STEP 1: importing the files. 
orig_GPR_file = pd.read_csv("RECON3_GPR_table_for_CORDA_full_V3.csv") #this is consistent. V3 is for typos.

#new addition: we have one Rx that is too large to have in a csv file, so it's in a text file. I'm going to try to add it manually.
text_file = open("Recon3_edited_ATPS4mi_Rx.txt", "r") #open text file in read mode
ATPS4mi_Rx = text_file.read() #read whole file to a string
text_file.close() #close file
del text_file

#substitute the string ATPS4mi_Rx for "HELPHELPHELP"
#print(orig_GPR_file.loc[orig_GPR_file['Rx_name'] == "ATPS4mi: adp[m] + 4.0 h[i] + pi[m] --> atp[m] + h2o[m] + 3.0 h[m]"])
orig_GPR_file= orig_GPR_file.replace("HELPHELPHELP", ATPS4mi_Rx)
#print(orig_GPR_file.loc[orig_GPR_file['Rx_name'] == "ATPS4mi: adp[m] + 4.0 h[i] + pi[m] --> atp[m] + h2o[m] + 3.0 h[m]"])

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
#keep only relevant columns
genename_lookup_table_abridged = genename_lookup_table_full[["name","geneName"]]
del genename_lookup_table_full

#second, import table for ccds_ids
ccds_table_full = pd.read_csv("CCDS_lookup_table.txt", delimiter = "\t")
ccds_table_abridged = ccds_table_full[["gene", 'gene_id', 'ccds_id']]
#the CCDS ids have "CCDS" in front of them; that's not helpful.
ccds_table_abridged["ccds_id"] = ccds_table_abridged["ccds_id"].str.replace('CCDS','')
del ccds_table_full

#Save the ccds table- we might use it later
#ccds_table_abridged.to_csv("gene_to_ccds_id_table.csv")
        
#STEP 2 converting RNAseq data to something useable
RNAseq_dict_ccds={}
RNAseq_dict_ccds_str={}

for key, value in RNAseq_dict.items():
    RNAseq_file_geneids = pd.merge(left=value, right= genename_lookup_table_abridged, how="left", left_on="target_id", right_on="name")
    RNAseq_file_geneids.drop(columns={"name"}, inplace=True)
    RNAseq_file_ccds = pd.merge(left=RNAseq_file_geneids, right=ccds_table_abridged, how="left", left_on="geneName", right_on="gene")
    RNAseq_file_ccds= RNAseq_file_ccds.drop(columns={"gene", "gene_id"})
    #RNAseq_dict_ccds.add(key, RNAseq_file_ccds)
    RNAseq_lookup_key = RNAseq_file_ccds.drop(columns={"target_id", "geneName"})
    RNAseq_lookup_key.drop_duplicates(subset= ["exp_score","ccds_id"], inplace=True)
    RNAseq_lookup_key = RNAseq_lookup_key.dropna(subset=['ccds_id'])
    RNAseq_lookup_key_nodup = RNAseq_lookup_key.sort_values('exp_score', ascending=False).drop_duplicates('ccds_id').sort_index()
    RNAseq_lookup_key_dict = RNAseq_lookup_key_nodup.set_index('ccds_id')['exp_score'].to_dict()
    keys_values = RNAseq_lookup_key_dict.items()
    RNAseq_lookup_key_dict_str = {str(key): str(value) for key, value in keys_values}
    #RNAseq_dict_ccds.update({key: RNAseq_file_ccds})
    RNAseq_dict_ccds.update({key: RNAseq_lookup_key_dict})
    RNAseq_dict_ccds_str.update({key: RNAseq_lookup_key_dict_str})
del key, value, RNAseq_file_geneids, RNAseq_file_ccds,RNAseq_lookup_key,RNAseq_lookup_key_nodup,RNAseq_lookup_key_dict,keys_values
del RNAseq_lookup_key_dict_str, RNAseq_dict_ccds
    
#word_freq.update({'before': 23}) 
    
#######STEP #3: doing the substitutions. 
#will need a dictionary of GPR_files and sample_names
def substituting_gprs(RNAseq_lookup_key_dict_str, sample_name):
    subrulecol = []
    sample_GPR_file = copy.deepcopy(orig_GPR_file)
    for row in sample_GPR_file["Rx_rule"]:
        if "(" in row:
            rowlist= re.split("[,*]", row) #start by comma and asterisk-separating items; ",|*" means "split by , or *
            rowlist2 = []
            for item in rowlist:
                try: 
                    float(item)
                    if item in RNAseq_lookup_key_dict_str.keys():
                        item2 = [val for key, val in RNAseq_lookup_key_dict_str.items() if item==key]
                        item3= item2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
                            #print(item2) 
                        rowlist2.append(item3)
                    else:
                        rowlist2.append("-1") #-1 means this ccds is not in the sample GPR file, so should not be included in the model
                except ValueError:
                    rowlist2.append(item) #this is where the items can't be converted to float (i.e. max/min)
            rowstring2 = ""
            rowstring2= ",".join([str(item3) for item3 in rowlist2]) #reassemble
            rowstring3 = rowstring2.replace("(,", "(")
            rowstring3= rowstring3.replace(",)",")")
            subrulecol.append(rowstring3)
        else:
            if row in RNAseq_lookup_key_dict_str.keys():
                 row2 = [val for key, val in RNAseq_lookup_key_dict_str.items() if row==key]
                 row3 = row2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
                 subrulecol.append(row3)
            else:
                subrulecol.append("-1")
    del row, rowlist, rowlist2, item, item2, item3, rowstring2, rowstring3, row2, row3              
    sample_GPR_file["edited_Rx_rule"]= subrulecol
    sample_GPR_file.to_csv("_GPRs_binned.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned.csv"):
        	dst=str(sample_name)+file
        	os.rename(src,dst)   
    return

#another option for this function- assigning 0s insteadof -1s
def substituting_gprs_round2(RNAseq_lookup_key_dict_str, sample_name):
    subrulecol = []
    sample_GPR_file = copy.deepcopy(orig_GPR_file)
    for row in sample_GPR_file["Rx_rule"]:
        if "(" in row:
            rowlist= re.split("[,*]", row) #start by comma and asterisk-separating items; ",|*" means "split by , or *
            rowlist2 = []
            for item in rowlist:
                try: 
                    float(item)
                    if item in RNAseq_lookup_key_dict_str.keys():
                        item2 = [val for key, val in RNAseq_lookup_key_dict_str.items() if item==key]
                        item3= item2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
                            #print(item2) 
                        rowlist2.append(item3)
                    else:
                        rowlist2.append("0") #try 0 instead of -1, as not enough Rxs being included
                except ValueError:
                    rowlist2.append(item) #this is where the items can't be converted to float (i.e. max/min)
            rowstring2 = ""
            rowstring2= ",".join([str(item3) for item3 in rowlist2]) #reassemble
            rowstring3 = rowstring2.replace("(,", "(")
            rowstring3= rowstring3.replace(",)",")")
            subrulecol.append(rowstring3)
        else:
            if row in RNAseq_lookup_key_dict_str.keys():
                 row2 = [val for key, val in RNAseq_lookup_key_dict_str.items() if row==key]
                 row3 = row2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
                 subrulecol.append(row3)
            else:
                subrulecol.append("0") #try 0 instead of -1, as not enough Rxs being included
    del row, rowlist, rowlist2, item, item2, item3, rowstring2, rowstring3, row2, row3              
    sample_GPR_file["edited_Rx_rule"]= subrulecol
    sample_GPR_file.to_csv("_GPRs_binned_moreperm.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned_moreperm.csv"):
        	dst=str(sample_name)+file
        	os.rename(src,dst)   
    return

        
######now, run it
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/")
for key, value in RNAseq_dict_ccds_str.items():
    substituting_gprs(value, key)

for key, value in RNAseq_dict_ccds_str.items():
    substituting_gprs_round2(value, key)

#######START HERE!!!

#DONE!!!!!!!!!!!
#okay, now CLEAR OUT all the files and re-import to see if we can get those equations solved.
binned_file_list = []
binned_file_names = []
for file in os.listdir():
    if("_binned_moreperm.csv" in file): #CHANGE AS NEEDED
        binned_file = pd.read_csv(file)
        binned_file.drop(columns={"Unnamed: 0", "Rx_rule"}, inplace=True)
        binned_file_list.append(binned_file)
        binned_file_names.append(str(file))
del binned_file, file

binned_file_names2 = [x.replace('_GPRs_binned_moreperm.csv', '') for x in binned_file_names] #CHANGE AS NEEDED
binned_dict = dict(zip(binned_file_names2, binned_file_list))
del binned_file_names

#import a dataframe of reactions without associated genes/GPRs
nongene_Rxs = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/RECON3_solved_GPRs_for_Rxs_without_genes_V2.csv")

#testing area
testo= binned_file_list[0]
Rx_solved=[]
for row in testo["edited_Rx_rule"]:
    try:
        solvedval = (eval(row)) #THE MAGIC TRICK! evaluates max/min statements. only works if max/min has been filled out correctly.
        Rx_solved.append(solvedval)
    except SyntaxError: #raised if max/min has a ccid still hanging out, or is missing a comma or something
        Rx_solved.append(row)
    except TypeError: #raised if there's already a single value (e/g 0) that's being read as a string
        Rx_solved.append(row)
testo["solved_Rx_rule"]= Rx_solved 
testo2 = pd.concat([testo, nongene_Rxs], axis=0)       #ah, nice- now we have 10600 Rxs, which is the total # in the Recon3 model

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
    value2.to_csv("_GPRs_binned_solved_V2.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned_solved_V2.csv"):
        	dst=str(key)+file
        	os.rename(src,dst)   

#great. let's move on.            

########################################Phfew. still running into issues. Let's try to look at this again.
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
test_file = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/s94_S90_GPRs_binned_solved.csv")
test_file.drop(columns={"Unnamed: 0"}, inplace=True)

ccds_table = pd.read_csv("gene_to_ccds_id_table.csv")
ccds_table.drop(columns={"Unnamed: 0"}, inplace=True)

#in the test_file, the edited_Rx_rule is the column we're interested in.
#where did edited_Rx_rule column come from? Look at the original file-edit there. 

    

##################################################EXTRA STUFF#######################################    
##################################################EXTRA STUFF#######################################    
##################################################EXTRA STUFF#######################################    


#conditionally adding commas to rowstring
        for i in rowlist2:
            rowstring2 = rowstring2 + str(i)
            if len(rowlist2)>=4 and i == rowlist2[0] or i == rowlist2[-1] or i == rowlist2[-2]:
                rowstring2 = rowstring2 + ""
            elif len(rowlist2)<4 and i == rowlist2[-1]:
                rowstring2 = rowstring2 + ""
            else:
                rowstring2 = rowstring2 + ","


for row in test_GPR_file1["Rx_rule"]:
    subrulecol = []
    if "(" in row:
        print(row)
        rowlist= re.split("[,*]", row) #start by comma and asterisk-separating items; ",|*" means "split by , or *
        #for item in rowlist:
            #item2 = re.split("*", item)
            #print(item2)
        #print(rowlist)
        rowlist2 = []
        for item in rowlist:
            if item in test_RNAseq_lookup_key_dict_str.keys():
                try:
                    float(item) #find out if item can be converted to float (i.e. is not max())
                    #print(item)
                    item2 = [val for key, val in test_RNAseq_lookup_key_dict_str.items() if item==key]
                    item3= item2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string
                    #print(item2) 
                    rowlist2.append(item3)
                except ValueError:
                    rowlist2.append(item)
            elif isinstance(item, float):
                rowlist2.append("0")
            else:
                rowlist2.append(item) #append string items (e.g. 'max(')
                #print(item)
        rowstring2 = ""
        print(rowlist2)
        rowstring2= ",".join([str(item3) for item3 in rowlist2]) #reassemble
        rowstring3 = rowstring2.replace("(,", "(")
        rowstring3= rowstring3.replace(",)",")")
        print(rowstring2)
        print(rowstring3)
        print("     ")
        subrulecol.append(rowstring3)
    else:
        if row in test_RNAseq_lookup_key_dict_str.keys():
             row2 = [val for key, val in test_RNAseq_lookup_key_dict_str.items() if row==key]
             row3= row2[0] #since the dictionary lookup returns a list (can't seem to make it work otherwise), convert to string 
