#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 09:56:43 2022

@author: adamo010
"""

#the goal here is to create a space where I can toggle around with reaction bins and such to get CORDA to work

#pulling code from both 01.20.22_combining_GPRs_and_expression_bins_V2.py and 01.26.22_testing_CORDA_with_Burns_data.py to test everythign at once. 

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
#substitute the string ATPS4mi_Rx for "HELPHELPHELP"
orig_GPR_file= orig_GPR_file.replace("HELPHELPHELP", ATPS4mi_Rx)
#import ONE test file
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto/")
RNAseq_file = pd.read_csv("s96_S92_expression_binned.csv")
RNAseq_file.drop(columns={"Unnamed: 0", "length", "eff_length", "est_counts", "tpm"}, inplace=True)
RNAseq_file_name = "s96_S92"
RNAseq_dict = {}
RNAseq_dict[RNAseq_file_name] = RNAseq_file #we're creating a dictionary here oNLY because it makes downstream stuff easier;
#we're just testing one file. See 01.20.22_combining_GPRs_and_expression_bins.py for full forloop
del RNAseq_file_name, RNAseq_file
#import table for matching geneids
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/")
genename_lookup_table_full = pd.read_csv("GRCh38_lookup_table.txt", delimiter = "\t")    
genename_lookup_table_abridged = genename_lookup_table_full[["name","geneName"]] #keep only relevant columns
del genename_lookup_table_full
#second, import table for ccds_ids
ccds_table_full = pd.read_csv("CCDS_lookup_table.txt", delimiter = "\t")
ccds_table_abridged = ccds_table_full[["gene", 'gene_id', 'ccds_id']]
ccds_table_abridged["ccds_id"] = ccds_table_abridged["ccds_id"].str.replace('CCDS','') #the CCDS ids have "CCDS" in front of them; that's not helpful.
del ccds_table_full

##############STEP 2 converting RNAseq data to something useable- basically, exchange weird gene ids into ccds
#RNAseq_dict_ccds={}
RNAseq_dict_ccds_str={} #I can't recall if we need both, but maybe?

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
    #RNAseq_dict_ccds.update({key: RNAseq_file_ccds})
    #RNAseq_dict_ccds.update({key: RNAseq_lookup_key_dict})
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
    sample_GPR_file.to_csv("_GPRs_binned_Feb2testo.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned_Feb2testo.csv"):
        	dst=str(sample_name)+file
        	os.rename(src,dst)   
    return
##now, run it
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/")
for key, value in RNAseq_dict_ccds_str.items():
    substituting_gprs(value, key)

###########Step 4: adding non-gene GPRs
binned_file_list = []
binned_file_names = []
for file in os.listdir():
    if("_binned_Feb2testo.csv" in file): #CHANGE AS NEEDED
        binned_file = pd.read_csv(file)
        binned_file.drop(columns={"Unnamed: 0", "Rx_rule"}, inplace=True)
        binned_file_list.append(binned_file)
        binned_file_names.append(str(file))
del binned_file, file

binned_file_names2 = [x.replace('_GPRs_binned_Feb2testo.csv', '') for x in binned_file_names] #CHANGE AS NEEDED
binned_dict = dict(zip(binned_file_names2, binned_file_list))
del binned_file_names
#import a dataframe of reactions without associated genes/GPRs. NOTE HERE: in V3, all of these have been assigned high confidence. The rationale
#is that these Rxs are part of the Recon3 model but have no associated genes, which means they must be essential for model running.
nongene_Rxs = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/RECON3_solved_GPRs_for_Rxs_without_genes_V3.csv")
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
    value2.to_csv("_GPRs_binned_solved_Feb2testo.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned_solved_Feb2testo.csv"):
        	dst=str(key)+file
        	os.rename(src,dst)   
del key, value, Rx_solved, row, solvedval, value2, file, dst, src

###################PART 2:::: parameterising the model####################
#now we are moving to code from 01.26.22_testing_CORDA_with_Burns_data.py
#use the magic command %reset in console to clear out environment.

import os
import corda
import copy
import pandas as pd
from corda import reaction_confidence
from cobra.io import load_matlab_model
from cobra.io import read_sbml_model
from corda import CORDA    

############step 5: import and clean up data. 
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models")
recon3 = load_matlab_model("/Users/adamo010/Documents/CORDA_CRC_human_models/Recon3D_301/Recon3DModel_301.mat")
#recon3 is the generic human model that we'll be parameterizing.

#import our test GPR tile
s96_S92_GPR_file= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/s96_S92_GPRs_binned_solved_Feb2testo.csv")
s96_S92_GPR_file.drop(columns={"Unnamed: 0", "edited_Rx_rule"}, inplace=True)
#also (maybe?) need to split Rx_name column
abb_Rx_name = []
for row in s96_S92_GPR_file["Rx_name"]:
    splitlist = row.split(r": ") #remove everything after :
    newrow = splitlist[0]
    abb_Rx_name.append(newrow)
del row, splitlist, newrow
s96_S92_GPR_file["abb_Rx_name"] = abb_Rx_name
s96_S92_GPR_file.drop(columns=['Rx_name'], inplace=True)     
gene_conf = dict(zip(s96_S92_GPR_file.abb_Rx_name, s96_S92_GPR_file.solved_Rx_rule))


############step 6: parameterize recon3 with RNAseq-based GPRs
#follow the tutorial here: https://github.com/resendislab/corda

#CORDA-ize the model
recon3_test1 = copy.deepcopy(recon3)
opt_test6 = CORDA(recon3_test1, gene_conf)
opt_test6.build() #began 11:02am. Ah, checked back at 2:24pm and it's done. 
#AND NOW WE WAIT

#once done, check # Rxs included etc.
print(opt_test6) #oKAY..., looking better- 6029 Rxs included

#see if model runs
opt_test6_model = opt_test6.model
opt_test6_model.objective="biomass_maintenance"
print(opt_test6_model.summary())
#well..... fuck. 

#what does RECON3 have going that test6 doesn't?
recon3.slim_optimize()
opt_test6_model.slim_optimize()

#what is the objective functon n recon3?
print(recon3.objective.expression)
print(opt_test6_model.objective.expression)

#same. Okay, let's try seeing which Rxs are present in recon3 but are missing in opt_test6 (based on a gapfilling thing)
from cobra.flux_analysis import gapfill

testing_gaps = gapfill(opt_test6_model, recon3, demand_reactions=False)
for reaction in testing_gaps[0]:
    print(reaction.id)

#well... "RuntimeError: failed to validate gapfilled model, try lowering the integer_threshold"
#let's try to validate the model.
import tempfile
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model
with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
    write_sbml_model(opt_test6_model, filename=f_sbml.name)
    report = validate_sbml_model(filename=f_sbml.name)
pprint(report)

#hmmm. "LibSBML error code -3: The requested action could not be performed. This can occur in a variety of contexts, such as 
#passing a null object as a parameter in a situation where it does not make sense to permit a null object.

#let's look at a few of our model components
print("exchanges", opt_test6_model.exchanges)
print("demands", opt_test6_model.demands)
print("sinks", opt_test6_model.sinks)
#loads of all of these. They're probaly the gene-free Rxs we set to 3 in this iteration of model creation.

#I guess we actually need to run gap filling.

#Hmm, lets just... try something real quick.
test_cobra_model_build = copy.deepcopy(opt_test6)
test_cobra_model_build.cobra_model()
print(test_cobra_model_build)
opt_test7_model = test_cobra_model_build.model
print(opt_test7_model.summary())
#ALLL RIGHT WHAT THE HELL?!?!?!?!
#while I am glad that we finally have some success, I am upset that it took this long

###############################################REEAL quick, I'm going to redo this with non-gene GPRs at 2 instead of 3. See if it still runs.
%reset

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
#substitute the string ATPS4mi_Rx for "HELPHELPHELP"
orig_GPR_file= orig_GPR_file.replace("HELPHELPHELP", ATPS4mi_Rx)
#import ONE test file
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto/")
RNAseq_file = pd.read_csv("s96_S92_expression_binned.csv")
RNAseq_file.drop(columns={"Unnamed: 0", "length", "eff_length", "est_counts", "tpm"}, inplace=True)
RNAseq_file_name = "s96_S92"
RNAseq_dict = {}
RNAseq_dict[RNAseq_file_name] = RNAseq_file #we're creating a dictionary here oNLY because it makes downstream stuff easier;
#we're just testing one file. See 01.20.22_combining_GPRs_and_expression_bins.py for full forloop
del RNAseq_file_name, RNAseq_file
#import table for matching geneids
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/")
genename_lookup_table_full = pd.read_csv("GRCh38_lookup_table.txt", delimiter = "\t")    
genename_lookup_table_abridged = genename_lookup_table_full[["name","geneName"]] #keep only relevant columns
del genename_lookup_table_full
#second, import table for ccds_ids
ccds_table_full = pd.read_csv("CCDS_lookup_table.txt", delimiter = "\t")
ccds_table_abridged = ccds_table_full[["gene", 'gene_id', 'ccds_id']]
ccds_table_abridged["ccds_id"] = ccds_table_abridged["ccds_id"].str.replace('CCDS','') #the CCDS ids have "CCDS" in front of them; that's not helpful.
del ccds_table_full

##############STEP 2 converting RNAseq data to something useable- basically, exchange weird gene ids into ccds
#RNAseq_dict_ccds={}
RNAseq_dict_ccds_str={} #I can't recall if we need both, but maybe?

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
    #RNAseq_dict_ccds.update({key: RNAseq_file_ccds})
    #RNAseq_dict_ccds.update({key: RNAseq_lookup_key_dict})
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
    sample_GPR_file.to_csv("_GPRs_binned_Feb2testo.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned_Feb2testo.csv"):
        	dst=str(sample_name)+file
        	os.rename(src,dst)   
    return
##now, run it
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/")
for key, value in RNAseq_dict_ccds_str.items():
    substituting_gprs(value, key)

###########Step 4: adding non-gene GPRs
binned_file_list = []
binned_file_names = []
for file in os.listdir():
    if("_binned_Feb2testo.csv" in file): #CHANGE AS NEEDED
        binned_file = pd.read_csv(file)
        binned_file.drop(columns={"Unnamed: 0", "Rx_rule"}, inplace=True)
        binned_file_list.append(binned_file)
        binned_file_names.append(str(file))
del binned_file, file

binned_file_names2 = [x.replace('_GPRs_binned_Feb2testo.csv', '') for x in binned_file_names] #CHANGE AS NEEDED
binned_dict = dict(zip(binned_file_names2, binned_file_list))
del binned_file_names
#import a dataframe of reactions without associated genes/GPRs. NOTE HERE: in V3, all of these have been assigned high confidence. The rationale
#is that these Rxs are part of the Recon3 model but have no associated genes, which means they must be essential for model running.
nongene_Rxs = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/RECON3_solved_GPRs_for_Rxs_without_genes_V3.csv") #EDITED HERE!!!
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
    value2.to_csv("_GPRs_binned_solved_Feb2testo.csv")
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_GPRs_binned_solved_Feb2testo.csv"):
        	dst=str(key)+file
        	os.rename(src,dst)   
del key, value, Rx_solved, row, solvedval, value2, file, dst, src

###################PART 2:::: parameterising the model####################
#now we are moving to code from 01.26.22_testing_CORDA_with_Burns_data.py
#use the magic command %reset in console to clear out environment.

import os
import corda
import copy
import pandas as pd
from corda import reaction_confidence
from cobra.io import load_matlab_model
from cobra.io import read_sbml_model
from corda import CORDA    

############step 5: import and clean up data. 
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models")
recon3 = load_matlab_model("/Users/adamo010/Documents/CORDA_CRC_human_models/Recon3D_301/Recon3DModel_301.mat")
#recon3 is the generic human model that we'll be parameterizing.

#import our test GPR tile
s96_S92_GPR_file= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/s96_S92_GPRs_binned_solved_Feb2testo.csv")
s96_S92_GPR_file.drop(columns={"Unnamed: 0", "edited_Rx_rule"}, inplace=True)
#also (maybe?) need to split Rx_name column
abb_Rx_name = []
for row in s96_S92_GPR_file["Rx_name"]:
    splitlist = row.split(r": ") #remove everything after :
    newrow = splitlist[0]
    abb_Rx_name.append(newrow)
del row, splitlist, newrow
s96_S92_GPR_file["abb_Rx_name"] = abb_Rx_name
s96_S92_GPR_file.drop(columns=['Rx_name'], inplace=True)     
gene_conf = dict(zip(s96_S92_GPR_file.abb_Rx_name, s96_S92_GPR_file.solved_Rx_rule))

############step 6: parameterize recon3 with RNAseq-based GPRs
#follow the tutorial here: https://github.com/resendislab/corda

#CORDA-ize the model
recon3_test1 = copy.deepcopy(recon3)
opt_test8 = CORDA(recon3_test1, gene_conf)
opt_test8.build() #began 7:51am
print(opt_test8) #how many Rxs included? 6029
opt_test8.cobra_model()
otp_test8_model = opt_test8.model
opt_test8.model.objective="biomass_maintenance" #not exactly sure where this goes, but hopefully somewhere
print(opt_test8.model.summary())

##########NEXT STEPS:: saving the models. 
cobra.io.save_matlab_model(opt_test8.model, "s96_S92_test_model.mat")

test_reload = cobra.io.load_matlab_model("s96_S92_test_model.mat")
print(test_reload.summary())
#nice.Works great. 
