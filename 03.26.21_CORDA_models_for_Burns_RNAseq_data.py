#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 13:58:47 2021

@author: adamo010
"""
#now that we've successfully (or at least attemptedly) binned the Burns RNAseq data into the appropriate categories, it's time to give
#building models the ol' college try. It's Friday at 2pm, so we'll see how far we get today. 

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

#let's start by importing the base model, which I've downloaded from vmh.life
recon3 = load_matlab_model("/Users/adamo010/Documents/CORDA_CRC_human_models/Recon3D_301/Recon3DModel_301.mat")
#recon3 is the generic human model that we'll be parameterizing.

len(recon3.reactions)
#this shows how many reactions are present. 10,600.

#To run CORDA, we need two things: a model and a set of confidences in the reactions.
#what I currently have is a set of confidences in the genes. So I'll need to convert that to confidences in reactions
#which will be... something. 


#let's start with a test file
#so, we have the gene_confidence rules; we also need the reaction rules

GPR_dict = cobra.manipulation.get_compiled_gene_reaction_rules(recon3)
#this gets us a dictionary of reactions and some _ast.Expression object that PRESUMABLY has the boolean rules in it. 

#okay, this code should extract and export the GPRs from recon3 to a csv
new_GPR_dict = {}

for key, value in GPR_dict.items():
    Rx= cobra.core.gene.ast2str(value)  #I'm used to get strings out of the ast values for the GPRs
    name= key
    new_GPR_dict.update({name:Rx})
    
with open('Recon3_GPR_dict.csv', 'w') as f:
    for key in new_GPR_dict.keys():
        f.write("%s,%s\n"%(key,new_GPR_dict[key]))

for key, value in new_GPR_dict.items():
    print(type(key))

del f, key, name, Rx, value

#Why the fuck is this so hard?
#try running a smaller model. 
import cobra.test
import os
from os.path import join

data_dir = cobra.test.data_dir

print("mini test files: ")
print(", ".join(i for i in os.listdir(data_dir) if i.startswith("mini")))

textbook_model = cobra.test.create_test_model("textbook")
ecoli_model = cobra.test.create_test_model("ecoli")
salmonella_model = cobra.test.create_test_model("salmonella")

ecoli_dict = cobra.manipulation.get_compiled_gene_reaction_rules(ecoli_model)

new_ecoli_dict = {}

for key, value in ecoli_dict.items():
    Rx= cobra.core.gene.ast2str(value)  #I'm used to get strings out of the ast values for the GPRs
    name= key
    new_ecoli_dict.update({name:Rx})

ecoli_gene = ecoli_model.genes.get_by_id("b4025")
ecoli_gene
print(ecoli_model.gene_reaction_rule)

all_EC_genes = ecoli_model.gene_reaction_rule

###########okay. let's give something a shot.

#1) change all dict values from and/or to max/min: dict is new_ecoli_dict
#2)run the following code:
#for key, value in new_ecoli_dict.items():
    #key.gene_reaction_rule = value
    #model.add_reactions

#step 1: going to export and re-import the new_ecoli_dict b/c I'm not sure how to edit strings in a dict and 
#honestly, i don't have the patience to fuck with it today. 

with open('Ecoli_test_GPR_dict.csv', 'w') as f:
    for key in new_ecoli_dict.keys():
        f.write("%s,%s\n"%(key,new_ecoli_dict[key]))
del f, key, name

#re-import
with open('Ecoli_test_GPR_dict_updated.csv', mode='r') as infile:
    reader = csv.reader(infile)
    with open('coors_new.csv', mode='w') as outfile:
        writer = csv.writer(outfile)
        updated_Ecoli_dict = {rows[0]:rows[1] for rows in reader}
del infile, reader, writer

#2)run the following code:
for key, value in updated_Ecoli_dict.items():
    value2= cobra.core.gene.parse_gpr(value)
    key.gene_reaction_rule = value2
    #print(key.gene_reaction_rule)
    ecoli_model.add_reactions

#get an AttributeError: 'str' object has no attribute 'gene_reaction_rule'
#I suppose this is not surprising, since there seems to be absolutely no way to change anything within the model itself.
#I think the issue is 'key' here. the 'key' replaces 'reaction', which is actually not a string.
#we need to figure out how to incorporate the model. Match the key to the reaction and replace the reaction with the value
#from the csv file.

#ideally, we'd iterate over all reactions in the model, and print them (as a start). But I do not know how to do that.
#let's try doing ONE reaction and converting it to max/min
#NI2uabcpp: atp_c + h2o_c + ni2_p --> adp_c + h_c + ni2_c + pi_c
# (b3479 and b3477 and b3476 and b3478 and b3480) or b4242
#should become max(min(b3479,b3477,b3476,b3478,b3480),b4242)

#first, print the rule.
NI2uabcpp = ecoli_model.reactions.get_by_id("NI2uabcpp") #this defines a reaction object, NI2uabcpp, from the ecoli_model
NI2uabcpp_gpr = NI2uabcpp.gene_reaction_rule #this extracts the gene reaction rule from the reaction object NI2uabcpp and saves it as a string
NI2uabcpp.genes #this prints a list of all genes involved in the reaction object NI2uabcpp_gpr 
NI2uabcpp.gene_reaction_rule = "max(min(b3479,b3477,b3476,b3478,b3480),b4242)" #define the new gene reaction rule
NI2uabcpp.gene_reaction_rule #prints the new gene reaction rule

#try something else
Rx_string_test = str("NI2uabcpp")
NI2uabcpp_p2 = ecoli_model.reactions.get_by_id(Rx_string_test) #this defines a reaction object, NI2uabcpp, from the ecoli_model


######building iterative scan:
for key, value in ecoli_dict.items():
    reaction = ecoli_model.reactions.get_by_id(str(key))
    reaction.gene_reaction_rule= str(value)

#okay, got some code. I guess now all that's left to do is hand-update the GPRs.

#check out some of the recon3 reactions that don't seem to have genes associated with them.
NDERSVitr =recon3.reactions.get_by_id("NDERSVitr")
NDERSVitr.gene_reaction_rule
PROFVShc= recon3.reactions.get_by_id("PROFVShc")
PROFVShc.gene_reaction_rule

#all right. I've made all the edits to the GPRs (switched from and/or to max/min.) Now I'd like to add these to the model

#first, import the updated rules
CORDA_Rx_rules_dict = pd.read_csv('Recon3_GPR_dict_CORDA_updated.csv', header=None, index_col=0, squeeze=True).to_dict()

#important note here: there are 10600 Rxs in GPR_dict, but only 2951 in CORDA_Rx_rules_dict; not every Rx in CORDA_Rx_rules_dict has
#an associated boolean. So, that might fuck things up a bit.

#well, try writing a foreloop to change this.
for key, value in CORDA_Rx_rules_dict.items():
    Rx_string = str(key)
    #print(Rx_string)
    rule_string = str(value)
    if key in new_GPR_dict:
        print(Rx_string)
        #Rx_name = recon3.reactions.get_by_id(Rx_string) #this defines a reaction object, Rx_name, from the recon3 model
        #name.gene_reaction_rule = rule_string #define the new gene reaction rule
   # else:
       # pass

for key, value in new_GPR_dict.items():
    print(key)
    print(value)

#running into an issue where the Reactions (keys) from my updated boolean csv file is not matching the reactions from the model.
for key, value in new_GPR_dict.items():
    print(type(key))

#aha, here might be thie issue. in new_GPR_dict.items, the keys (Rx names) aren't strings, they're cobra.core.reaction.Reaction objects
#convert to strings

string_GPR_dict = {}
for key, value in new_GPR_dict.items():
    key2= str(key)
    value2= str(value)
    string_GPR_dict.update({key2:value2})
    
for key, value in string_GPR_dict.items():
    print(key)
#better.

for key, value in CORDA_Rx_rules_dict.items():
    Rx_string = str(key)
    #print(Rx_string)
    rule_string = str(value)
    if key in string_GPR_dict:
        #print(type(Rx_string))
        #print(type(key))
        #if key is Rx_string:
            #print("FINE") #these fucking things should match
        Rx_name = recon3.reactions.get_by_id(Rx_string) #this defines a reaction object, Rx_name, from the recon3 model
        print(Rx_name)
        Rx_name.gene_reaction_rule = rule_string #define the new gene reaction rule
    else:
        pass

#try changing ONE FUCKING REACTION AT ALL in the recon3 model.
#<Reaction GALT at 0x7ff58aaadf90>: '2592.1 or 2592.3 or 2592.2'
GALT_Rx = recon3.reactions.get_by_id("GALT") 
GALT_Rx.gene_reaction_rule

for key, value in CORDA_Rx_rules_dict.items():
    Rx_string = str(key)
    print(Rx_string)
    rule_string = str(value)
    if key in string_GPR_dict:
        print(type(Rx_string))
        #print(type(key))
        if key is Rx_string:
            print("FINE") #these fucking things should match
        Rx_name = recon3.reactions.get_by_id(Rx_string) #this defines a reaction object, Rx_name, from the recon3 model
        print(Rx_name)
        Rx_name.gene_reaction_rule = rule_string #define the new gene reaction rule
    else:
        pass

#aha, I know what's happening. I need to separate CORDA_Rx_rules_dict by colons; 
#NADH2_u10mi: 5.0 h[m] + nadh[m] + q10[m] --> 4.0 h[i] + nad[m] + q10h2[m] isn't actually the Rx name, it's the Rx itself. 

#so, re-import the updated boolean rules file, trim off the Rxs from the first column, and re-make the dictionary. All else should be the same.
#first, import the updated rules
CORDA_Rx_rules2 = pd.read_csv('Recon3_GPR_dict_CORDA_updated.csv', header=None, index_col=0, squeeze=True).to_dict()
sep = ":"
CORDA_Rx_rules2 = {k.split(sep, 1)[0]: v for (k, v) in CORDA_Rx_rules2.items()} #nice. this removes all text afer the sep (:) value in the dictionary keys
del sep

#try again.
for key, value in CORDA_Rx_rules2.items():
    Rx_string = str(key)
    #print(Rx_string)
    new_rule_string = str(value)
    #print(new_rule_string)
    if key in string_GPR_dict:
        #print(type(Rx_string))
        #print(type(key))
        #if key is Rx_string:
            #print("FINE") #these fucking things should match
        Rx_name = recon3.reactions.get_by_id(Rx_string) #this defines a reaction object, Rx_name, from the recon3 model
        print(Rx_name)
        Rx_name.gene_reaction_rule = new_rule_string #define the new gene reaction rule
        Rx_name.gene_reaction_rule
    else:
        pass

#suspiciously fast. 
GALT_Rx = recon3.reactions.get_by_id("GALT") 
GALT_Rx.gene_reaction_rule
#b/c it didn't work

#the get_by_id step is where the issue is.
#PLEASE tell me I don't have to do all these by hand.
for key, value in CORDA_Rx_rules2.items():
    Rx_string = str(key)
    new_rule_string = str(value)
    if new_rule_string != "nan":
        #print(key, value)
        Rx_name = recon3.reactions.get_by_id(Rx_string) #this defines a reaction object, Rx_name, from the recon3 model
        Rx_name.genes #doesn't work
        Rx_name.gene_reaction_rule = new_rule_string #define the new gene reaction rule
        Rx_name.gene_reaction_rule

#
#FAOXC101C8x max(1576.1,max(4051.1,8529.1))
FAOXC101C8x_Rx = recon3.reactions.get_by_id("FAOXC101C8x")
FAOXC101C8x_Rx.genes
FAOXC101C8x_Rx.gene_reaction_rule
#well, it... seems like those are changed? did it actually work?

#here's the big question- does the model run at all?




















