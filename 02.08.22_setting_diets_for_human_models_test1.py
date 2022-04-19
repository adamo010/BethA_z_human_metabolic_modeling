#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 14:03:14 2022

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

#our goal here is to make sure our cool models run, and also that we can change the growth media

os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/")

#part 1: import test files
test_model = cobra.io.load_matlab_model("s96_S92_test_model.mat")

highfiber_diet = pd.read_csv("/Users/adamo010/Documents/microbiome_on_diets/highfiber_fluxes.tsv", sep= "\t")

#part 2: poke around in model. 
print(test_model.medium) #wow, okay, looks like the default is EVERYTHING IN EXCESS!
print(test_model.optimize()) #our default solutuion is 193.697. this is MAX growth

max_growth = test_model.slim_optimize() #save the max/default solution
minimal_medium(test_model, max_growth) #look at minimal medium; it's actually pretty small.

#well, let's give  it a shot.
#test_model.medium = highfiber_diet #nope.

#oh right, a dictionary works better here.
highfiber = {} #editing with every iteration, so need to reset.
with open("/Users/adamo010/Documents/microbiome_on_diets/highfiber_fluxes.tsv", "r") as file:
    next(file)
    for line in file:                                   
        linestuff = line.strip().split()                    
        highfiber[linestuff[0]] = float(linestuff[1])  #this creates a dictionary called 'highfiber' containing the metabolite (key) and its amount (value) from the supplied diet file
highfiber_edited = {k.replace("[e]", "_m"): v for k, v in highfiber.items()}
highfiber_edited2 = {k.replace("(e)", "_m"): v for k, v in highfiber_edited.items()}
del file, line, linestuff

test_model.medium= highfiber #used to be master_MM, but now I'm generating locally

#FUCK's sake, fine, I'll go through the many hoops.
#pull model MM
model_only_MM =  minimal_medium(test_model, max_growth)
model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's coop tradeoff minimal Rxs and their fluxes

#then, find shared metabolites between diet and metabolites that model can import
test_model_exchange_ids = [exchange.id for exchange in test_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
#identify metabolites from EUstandard diet that can be imported
to_keep_from_diet = {}
for key, value in highfiber.items():
    if key in test_model_exchange_ids:
    	to_keep_from_diet[key]= value
del key, value
###merge these dicts
keys = set(to_keep_from_diet.keys()).union(model_only_MM_dict.keys())
sample_driven_diet_amounts = {k:max(to_keep_from_diet.get(k,float('-inf')), model_only_MM_dict.get(k, float('-inf'))) for k in keys}
#set model MM
test_model.medium= sample_driven_diet_amounts #ding   
sample_sol = test_model.optimize()
fluxes = sample_sol.fluxes
test_model_summary = print(test_model.summary()) #this is so nice. How do I... save this?

with open('filename.txt', 'w') as f:
    print(test_model.summary(), file=f) #this saves it as a text file. useful, but I'd really like a csv file. 
text_file = open("filename.txt", "r")
model_summary_text = text_file.read() #read whole file to a string
text_file.close() #close file

model_summary_list = model_summary_text.split("------") #close, but not quite
model_summary_list2 = model_summary_text.splitlines() #now every line is a list element.
model_summary_list3 = model_summary_text.split("\n\n") #that's what I want- split by empty line

#okay, now I should be able to dump each element in this list into a csv: element 1 is Objective/biomass,
#element 2 is Update, element 3 is secretion

#probaby need to do them all individually.
text_file = open("objective.txt", "w")
n = text_file.write(model_summary_list3[0])
text_file.close()
objective_df = pd.read_csv("objective.txt")
#that one is fine, I guess.

#try a more complicated one.
text_file = open("uptake.txt", "w")
n = text_file.write(model_summary_list3[1])
text_file.close()
uptake_df = pd.read_csv("uptake.txt", sep= "\t") #START HERE IT WONT READ AS TABLE
#first, remove first row; it's just ----- spacers anyway
uptake_df = uptake_df.iloc[1: , :]
#now, split into five columns
uptake_df[['Metabolite','Reaction', 'Flux', "C-number", "C-Flux"]] = uptake_df.Uptake.str.split(expand=True)
#okay, that... mostly, worked, I guess. now delete first column (which is the old string column) and first row (which is now column headers)
uptake_df = uptake_df.iloc[1: , :] #delete first row
uptake_df = uptake_df.iloc[: , 1:] #delete first column
uptake_df.to_csv("test_sample_uptake_df.csv")

#real quick, let's make sure secretion works too.
text_file = open("secretion.txt", "w") #note that secretion.txt is a placeholder, and will be deleted. 
text_file.write(model_summary_list3[2])
text_file.close()
secretion_df = pd.read_csv("secretion.txt", sep= "\t") #START HERE IT WONT READ AS TABLE
#first, remove first row; it's just ----- spacers anyway
secretion_df = secretion_df.iloc[1: , :]
    #now, split into five columns
secretion_df[['Metabolite','Reaction', 'Flux', "C-number", "C-Flux"]] = secretion_df.Secretion.str.split(expand=True)
    #okay, that... mostly, worked, I guess. now delete first column (which is the old string column) and first row (which is now column headers)
secretion_df = secretion_df.iloc[1: , :] #delete first row
secretion_df = secretion_df.iloc[: , 1:] #delete first column
    secretion_df.to_csv("_V01_Burns2015_human_highfiber_uptake.csv")
    os.remove("secretion.txt")






    
    
#well... cool. I guess let's do this in a big function?
#I think a function rather than a foreloop is the way to go, because we can specify both sample file and sample id.
def running_human_models(model_file, sample_id):
    os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_sample_specific_human_models/")
    test_model = cobra.io.load_matlab_model(model_file) #import the matlab model file
    model_only_MM =  minimal_medium(test_model, max_growth) #get the model's MM
    model_only_MM_dict = model_only_MM.to_dict() #creates a dictionary of model's minimal Rxs and their fluxes
    os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/highfiber_diet_model_outputs/")
    highfiber = {} #initiate a diet dictionary of metabolites and their amounts
    with open("/Users/adamo010/Documents/microbiome_on_diets/highfiber_fluxes.tsv", "r") as file:
        next(file)
        for line in file:                                   
            linestuff = line.strip().split()                    
            highfiber[linestuff[0]] = float(linestuff[1])  #this creates a dictionary called 'highfiber' containing the metabolite (key) and its amount (value) from the supplied diet file
    del file, line, linestuff
    test_model_exchange_ids = [exchange.id for exchange in test_model.exchanges]  #list comprehension: gets all the export Rxs (exchange ids) from the  model   
    to_keep_from_diet = {} #initiate a dictionary of metabs/amounts to keep from highfiber dict
    for key, value in highfiber.items():
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
    fluxes.to_csv("_V01_Burns2015_human_highfiber_fluxes.csv")
    #step 2: save model summary. BUckle up, this is messy.
    with open('filename.txt', 'w') as f:
        print(test_model.summary(), file=f) #this saves it as a text file. useful, but I'd really like a csv file. 
    text_file = open("filename.txt", "r")
    model_summary_text = text_file.read() #read whole file to a string
    text_file.close() #close file
    model_summary_list3 = model_summary_text.split("\n\n") #that's what I want- split by empty line, into the three "tables" of useful info
    #####okay, now save each chunk of model_summary into three files: objective/biomass, Uptake, and Secretion
    ###start with the easy one: objective. Leave it as a text file; can edit that mess into something useful later
    text_file = open("_V01_Burns2015_human_highfiber_objective.txt", "w")
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
    uptake_df.to_csv("_V01_Burns2015_human_highfiber_uptake.csv")
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
    secretion_df.to_csv("_V01_Burns2015_human_highfiber_secretion.csv")
    os.remove("secretion.txt")
    #THE GREAT RENAMENING
    for file in os.listdir():
        src=file
        if fnmatch.fnmatch(file, "_V01_Burns2015_human_highfiber_fluxes.csv"):
            dst=str(sample_id)+file
            os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V01_Burns2015_human_highfiber_objective.txt"):
            dst=str(sample_id)+file
            os.rename(src,dst)
        elif fnmatch.fnmatch(file, "_V01_Burns2015_human_highfiber_uptake.csv"):
            dst=str(sample_id)+file
            os.rename(src,dst)    
        elif fnmatch.fnmatch(file, "_V01_Burns2015_human_highfiber_secretion.csv"):
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
    running_human_models  (value, key)


#well... it looks like not all the models made something useful.
#let's first try to download them all and at least get output.


##########have moved! to 02.09.22 version of this file












