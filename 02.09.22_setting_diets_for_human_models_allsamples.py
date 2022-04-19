#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 11:46:33 2022

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
from cobra.io import load_matlab_model
import cobra
from corda import CORDA
from cobra.medium import minimal_medium

#building off 02.08.22_setting_diets_for_human_models_test1.py

os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/")

#I think a function rather than a foreloop is the way to go, because we can specify both sample file and sample id.
def running_human_models(model_file, sample_id):
    os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/")
    test_model = cobra.io.load_matlab_model(model_file) #import the matlab model file
    max_growth = test_model.slim_optimize() #save the max/default solution; need it for later
    model_only_MM =  minimal_medium(test_model, max_growth) #get the model's MM
    model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
    os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/EUstd_diet_model_outputs/")
    EUstd = {} #initiate a diet dictionary of metabolites and their amounts
    with open("/Users/adamo010/Documents/microbiome_on_diets/EUstandard_fluxes.tsv", "r") as file:
        next(file)
        for line in file:                                   
            linestuff = line.strip().split()                    
            EUstd[linestuff[0]] = float(linestuff[1])  #this creates a dictionary called 'highfiber' containing the metabolite (key) and its amount (value) from the supplied diet file
    del file, line, linestuff
    test_model_exchange_ids = [exchange.id for exchange in test_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
    to_keep_from_diet = {} #initiate a dictionary of metabs/amounts to keep from highfiber dict
    for key, value in EUstd.items():
        if key in test_model_exchange_ids:
            to_keep_from_diet[key]= value
    del key, value
    keys = set(to_keep_from_diet.keys()).union(model_only_MM_dict.keys()) #identify metabolites from diet that can be imported by model
    sample_driven_diet_amounts = {k:max(to_keep_from_diet.get(k,float('-inf')), model_only_MM_dict.get(k, float('-inf'))) for k in keys} #merge these dicts
    #now, we're ready to set model medium and run the model
    test_model.medium= sample_driven_diet_amounts #ding   
    sample_sol = test_model.optimize() #rate limiting step; might take a while.
    #step 1: saving fluxes
    fluxes = sample_sol.fluxes
    fluxes.to_csv("_V01_Burns2015_human_EUstd_fluxes.csv")
    #step 2: save model summary. BUckle up, this is messy.
    with open('filename.txt', 'w') as f:
        print(test_model.summary(), file=f) #this saves it as a text file. useful, but I'd really like a csv file. 
    text_file = open("filename.txt", "r")
    model_summary_text = text_file.read() #read whole file to a string
    text_file.close() #close file
    model_summary_list3 = model_summary_text.split("\n\n") #that's what I want- split by empty line, into the three "tables" of useful info
    #####okay, now save each chunk of model_summary into three files: objective/biomass, Uptake, and Secretion
    ###start with the easy one: objective. Leave it as a text file; can edit that mess into something useful later
    text_file = open("_V01_Burns2015_human_EUstd_objective.txt", "w")
    text_file.write(model_summary_list3[0])
    text_file.close()
    del text_file
    ###next, something harder: Uptake
    text_file = open("uptake.txt", "w") #note that uptake.txt is a placeholder, and will be deleted. 
    text_file.write(model_summary_list3[1])
    text_file.close()
    uptake_df = pd.read_csv("uptake.txt", sep= "\t")
    #first, remove first row; it's just ----- spacers anyway
    uptake_df = uptake_df.iloc[1: , :]
    #now, split into five columns
    uptake_df[['Metabolite','Reaction', 'Flux', "C-number", "C-Flux"]] = uptake_df.Uptake.str.split(expand=True)
    #okay, that... mostly, worked, I guess. now delete first column (which is the old string column) and first row (which is now column headers)
    uptake_df = uptake_df.iloc[1: , :] #delete first row
    uptake_df = uptake_df.iloc[: , 1:] #delete first column
    uptake_df.to_csv("_V01_Burns2015_human_EUstd_uptake.csv")
    os.remove("uptake.txt")
    #finally, let's move on to Secretion
    text_file = open("secretion.txt", "w") #note that secretion.txt is a placeholder, and will be deleted. 
    text_file.write(model_summary_list3[2])
    text_file.close()
    secretion_df = pd.read_csv("secretion.txt", sep= "\t")
    #first, remove first row; it's just ----- spacers anyway
    secretion_df = secretion_df.iloc[1: , :]
    #now, split into five columns
    secretion_df[['Metabolite','Reaction', 'Flux', "C-number", "C-Flux"]] = secretion_df.Secretion.str.split(expand=True)
    #okay, that... mostly, worked, I guess. now delete first column (which is the old string column) and first row (which is now column headers)
    secretion_df = secretion_df.iloc[1: , :] #delete first row
    secretion_df = secretion_df.iloc[: , 1:] #delete first column
    secretion_df.to_csv("_V01_Burns2015_human_EUstd_secretion.csv")
    os.remove("secretion.txt")
    #THE GREAT RENAMENING
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_V01_Burns2015_human_EUstd_fluxes.csv"):
            dst=str(sample_id)+file
            os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V01_Burns2015_human_EUstd_objective.txt"):
            dst=str(sample_id)+file
            os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V01_Burns2015_human_EUstd_uptake.csv"):
            dst=str(sample_id)+file
            os.rename(src,dst)    
        elif fnmatch.fnmatch(file, "_V01_Burns2015_human_EUstd_secretion.csv"):
            dst=str(sample_id)+file
            os.rename(src,dst)    
    del test_model
    return
            

#now, let's make a dictionary of all these sample ids and their file names, and use that to run our function
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models")
human_model_list = []
human_model_names =[]
for file in os.listdir():
    if("_V3_model.mat" in file):
        human_model_list.append(file)
        human_model_names.append(str(file))
del file
human_model_names2 = [x.replace('_V3_model.mat', '') for x in human_model_names] #CHANGE AS NEEDED
model_file_dict = dict(zip(human_model_names2, human_model_list))
del human_model_names
  
#we'll start with a high fiber diet, just for shiggles.
for key, value in model_file_dict.items():
    running_human_models(value, key)

#well... it looks like not all the models made something useful.
#let's first try to download them all and at least get output.
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/")

for key, value in model_file_dict.items():
    test_model = cobra.io.load_matlab_model(value) #import the matlab model file
    max_growth = test_model.slim_optimize() #save the max/default solution; need it for later
    print(max_growth)
    #model_only_MM =  minimal_medium(test_model, max_growth) #get the model's MM
    #model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes

#OH HANG ON, did I miss this? max_growth = test_model.slim_optimize() #save the max/default solution
#add that to the function and retry. 
#now I have a ValueError- did not read any bytes. Me neither, program. Me neither. 

for key, value in model_file_dict.items():
    test_model = cobra.io.load_matlab_model(value) #import the matlab model file
    max_growth = test_model.slim_optimize() #save the max/default solution; need it for later
    print(max_growth)

#all right, it looks like we have a model that doesn't work. THe first 19 models in the dictionary printed max_growth.
#therefore, the issue must be with #20. But, since dictionaries are unordered, that won't be enough.
#good thing we have a dictionary and can print other things. like sample names
for key, value in model_file_dict.items():
    print(key)
    test_model = cobra.io.load_matlab_model(value) #import the matlab model file
    max_growth = test_model.slim_optimize() #save the max/default solution; need it for later
    print(max_growth)

#looks like it's s37_S59
ts_model = cobra.io.load_matlab_model("s37_S59_V3_model.mat")
#interesting... I wonder what that's about. That wasn't even one of the funky models that didn't generate the first time.

#let's try a different method of importing data
from os.path import dirname, join as pjoin
import scipy.io as sio
ts_model = sio.loadmat("s37_S59_V3_model.mat")

ts_model = sio.loadmat("B04_S4_V3_model.mat") #okay, so this one works. Clearly it's not the importing. 
#well, maybe something happened with the file creation or something. let's figure out which files are causing the problem and rebuild their models.

problem_sample_ids = []
for key, value in model_file_dict.items():
    #print(key)
    try:
        test_model = cobra.io.load_matlab_model(value) #import the matlab model file
    except ValueError:
        problem_sample_ids.append(key)
        
#Conveniently, there are only two: s37_S59, and B08_S8. Let's remake those.     
#remade! So, let's try them.
ts_model_1 = cobra.io.load_matlab_model("s37_S59_V3_model.mat") #okay great! It's importing!
ts_model_2 = cobra.io.load_matlab_model("B08_S8_V3_model.mat")

max_growth1 = ts_model_1.slim_optimize() #save the max/default solution; need it for later
max_growth2 = ts_model_2.slim_optimize() #save the max/default solution; need it for later
#thumbs up. Let's run everything again. 

        
