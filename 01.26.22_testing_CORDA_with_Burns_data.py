#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 15:00:04 2022

@author: adamo010
"""
#borrowing from 03.09.21_CORDA_tutorial.py; Thanks to Osbaldo Resendis for writing a python wrapper
#for a matlab version of CORDA: Cost Optimization Reaction Dependency Assessment. 
#This is from Amina Qutub's lab, and is many of the options for creating tissue-specific metabolic models.
#I am using it because a) the paper made something that is most similar to the kinds of models
#I want to generate, and b) b/c I saw Dr. Qutub speak at a conference and she seemed nice.

#we're working in CORDA_human_fba for our virtual environment
import os
import corda
import pandas as pd

#first, switch to our NEW FOLDER
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models")

from corda import reaction_confidence
from cobra.io import load_matlab_model

#import our corda model
recon3 = load_matlab_model("/Users/adamo010/Documents/CORDA_CRC_human_models/Recon3D_301/Recon3DModel_301.mat")
#recon3 is the generic human model that we'll be parameterizing.

len(recon3.reactions) #gives the number of reactions in the model; it is 10600.

#import a test GPR file
S94_s90_GPR_file= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/s94_S90_GPRs_binned_solved.csv")
S94_s90_GPR_file.drop(columns={"Unnamed: 0", "edited_Rx_rule"}, inplace=True)
#also (maybe?) need to split Rx_name column
abb_Rx_name = []
for row in S94_s90_GPR_file["Rx_name"]:
    splitlist = row.split(r": ") #remove everything after :
    newrow = splitlist[0]
    abb_Rx_name.append(newrow)
del row, splitlist, newrow
S94_s90_GPR_file["abb_Rx_name"] = abb_Rx_name

S94_s90_GPR_file.drop(columns=['Rx_name'], inplace=True)     

######################other stuff, nothing to see here.
#pull out one Rx
ATPS4mi_Rx = recon3.reactions.get_by_id("ATPS4mi") 
ATPS4mi_Rx.gene_reaction_rule

##########ok, convert GPR file to dictionary
gene_conf = dict(zip(S94_s90_GPR_file.abb_Rx_name, S94_s90_GPR_file.solved_Rx_rule))

#follow the tutorial here: https://github.com/resendislab/corda
from cobra.io import read_sbml_model
from corda import CORDA    

opt = CORDA(recon3, gene_conf)
#ooo! did it work?
print(opt)
print([opt.model.reactions.get_by_id(k).reaction for k, used in opt.included.items() if used])

#so.... how do I add the parameters?
opt = CORDA(recon3, gene_conf)

#according to CORDA_running_instructions.txt (more detail), those parameters are implicit. THe full command is as follows:
#CORDA(model,metTests,ES,PR,NP,PRtoNP,constraint,constrainby,om,ntimes,nl)
    #model: recon3 (always)
    #metTests/ES/PR/NP: covered in gene_conf (unique to each sample)- metabolic Rxs and their confidences
    #PRtoNP: threshold to include NP reactions. Default value of 2
    #constraint: constraint value. Default vaule of 1
    #constrainby: type of constraint used when defining reaction dependency. Two options: perc or val. Default val.
    #om: cost assigned to reactions while determining reaction dependency. Default value 1e+04
    #ntimes: # times simulations are performed to determine Rx dependecy. Default value of 5
    #nl: Noise level. Numeric. Sampled uniformly between zero and nl. Max default is 10-2

#long story short, I'm just going to have those default parameters and leave as is. 

#opt (the corda object) has a lot of features, only one of which is the model. Try to pull out the model.
s94_s90_model = opt.model

#saving models: historically I've used a pickle file, but the cobrapy docs recommend against that. Will try SMBL.
import cobra
cobra.io.write_sbml_model(opt.model, "S94_s90_human_model.xml")

#did not work. Maybe saving as a matlab model, since Recon3 is a matlab model?
cobra.io.save_matlab_model(s94_s90_model, "S94_s90_human_model.mat")
#okay. Re-import:
S94_s90_reimport_test = cobra.io.load_matlab_model("S94_s90_human_model.mat")

#seems to work. But does it run?
S94_s90_solution = S94_s90_reimport_test.optimize() #SUSPICIOUSLY fast
print(S94_s90_solution)
#yea, okay. GR=0

#let's try to troubleshoot this.
for x in S94_s90_reimport_test.reactions:
    print("%s : %s" % (x.id, x.reaction)) #looks fine.
for x in S94_s90_reimport_test.metabolites:
    print("%s : %s" % (x.id, x.formula)) #that's a list of metabolites and their chemical compositions!
    
print(S94_s90_reimport_test.objective.expression) #what the hell, does this model not have an objective function??????

#oh, look- here's some code for validating the model.
import tempfile
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model
with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
    write_sbml_model(s94_s90_model, filename=f_sbml.name)
    report = validate_sbml_model(filename=f_sbml.name)
pprint(report)

#whoops. Sure enough, no objective function. Is this a me issue or a recon issue?
with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
    write_sbml_model(recon3, filename=f_sbml.name)
    report = validate_sbml_model(filename=f_sbml.name)
pprint(report)
#recon3 doesn't throw up any errors, so it's probably a me issue. What's the objective in recon3?
print(recon3.objective.expression)
#1.0*biomass_maintenance - 1.0*biomass_maintenance_reverse_95d2f

#so... can we just set model_objective to biomass_maintenance?
s94_s90_model.objective="biomass_maintenance"
#seems to work. But does it run?
s94_s90_solution = s94_s90_model.optimize()
print(s94_s90_solution)
#nope.
with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
    write_sbml_model(s94_s90_model, filename=f_sbml.name)
    report = validate_sbml_model(filename=f_sbml.name)
pprint(report)
#looks better (i.e. no more errors), but still getting a solution of 0.00

#does the same thing happen with Recon3?
recon3_solution = recon3.optimize()
print(recon3_solution)
#NO! WHAT THE HELL?

print(s94_s90_model.summary())
print(recon3.summary())
#yeah..... broken.

#let's try rebuilding the model.
opt = CORDA(recon3, gene_conf)
opt.build() #somehow missed adding this in before- is it going to be the magic thing that fixes everything?
print(opt)
#EEEEE.... okay. THat seems to have fixed things. It did take 1000 years/20min to run, though.

s94_s90_model = opt.model #save model component
with tempfile.NamedTemporaryFile(suffix='.xml') as f_sbml:
    write_sbml_model(s94_s90_model, filename=f_sbml.name)
    report = validate_sbml_model(filename=f_sbml.name)
pprint(report)
#okay, still have to specify that biomass_maintenance is the objective; I guess that's okay. 
s94_s90_model.objective="biomass_maintenance"
#yep, will still need to save as a .mat file vs an smbl file, but so far, no errors.
s94_s90_solution = s94_s90_model.optimize()
print(s94_s90_solution)
#HUHHHHHHH
print(s94_s90_model.summary())
#still getting an empty dataframe.
#maybe try with opt, instead of model.
s94_s90_solution2 = opt.optimize() #nope, AttributeError: 'CORDA' object has no attribute 'optimize'

#I wonder if there's something wrong with the gene_conf file.
#hmmmm... the corda github from the resendis lab says that the geneconf file should be formatted as follows:
#gene_conf = {"gene1": 1, "gene2": 3, "gene4": -1} # missing entries are automatically assigned zeroes
#I.... made something like {"Rx_1": 1, "Rx_2": -1} # you know, that code that took two weeks to figure out?

#fuck it. Let's try building the model with the gene-bin values vs the RX-bin values.
#borrow from 01.20.22_combining_GPRs_and_expression_bins_V2.py
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
#ccds_table_abridged["ccds_id"] = ccds_table_abridged["ccds_id"].str.replace('CCDS','')
#del ccds_table_full
#STEP 2 converting RNAseq data to something useable
#the larger fcn is in the abovementioned file; I"m just doing one sample.

s94_s90_binned_RNAseq = pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto/s94_S90_expression_binned.csv")
s94_s90_binned_RNAseq.drop(columns={"Unnamed: 0", "length", "eff_length", "est_counts", "tpm"}, inplace=True)

s94_s90_file_geneids = pd.merge(left=s94_s90_binned_RNAseq, right= genename_lookup_table_abridged, how="left", left_on="target_id", right_on="name")
s94_s90_file_geneids.drop(columns={"name"}, inplace=True)
#s94_s90_file_ccds = pd.merge(left=s94_s90_file_geneids, right=ccds_table_abridged, how="left", left_on="geneName", right_on="gene")
#s94_s90_file_ccds= s94_s90_file_ccds.drop(columns={"gene", "gene_id"})
s94_s90_lookup_key = s94_s90_file_geneids.drop(columns={"target_id"})
s94_s90_lookup_key_nodup = s94_s90_lookup_key.sort_values('exp_score', ascending=False).drop_duplicates('geneName').sort_index()
s94_s90_lookup_key_dict = s94_s90_lookup_key_nodup.set_index('geneName')['exp_score'].to_dict()

#okay, try AGAIN.

opt_round3 = CORDA(recon3, s94_s90_lookup_key_dict)
#nope, that didn't work either. get that ValueError: 10FTHF5GLUtl missing from confidences! thing
#which means we're going to run into issues with reactions like that one, which have no associated genes. When we ran into this with the
#reaction binning, we could just bin the reaction. Now, with gene binning, that's not an option. Which means we need to figure out a way to bin
#the reactions instead of the genes. 

#okay, s94_s90_lookup_key_dict is the genes, and gene_conf is the reactions

for r in recon3.reactions:
    print(r)
    #print(conf[r.id])

#testing dummy model
len(recon3.reactions) #have 10600 reactions
conf = {}
for r in recon3.reactions: conf[r.id] = -1
conf["biomass_maintenance"] = 3
#set biomass as the only reaction that HAS to be included, all other Rxs as "not include" (-1)
#okay, so key is Rx, value is bin value.

opt_test = CORDA(recon3, conf)
opt_test.build()
print(opt_test)
opt_test_model = opt_test.model
opt_test_model.objective="biomass_maintenance"
print(opt_test_model.summary())
#I don't know if it's a good or a bad thing that this didn't work either. 

#while we wait for build to run, I am wondering waht the biomass Rx is in our s94_s90_lookup_key_dict. It should be 3; what if it's not?
for key, value in s94_s90_lookup_key_dict.items():
    if key == "biomass_maintenance":
        print(key, '::', value)
#hmmm. biomass_maintenance isn't there? Need to find out what our high confidence Rxs are, and add biomass to them if they're not there. 

#let's check the default Rx values for recon3 (remember to re-import, if you've edited)
for r in recon3.reactions:     
    print(conf[r.id])
#that's kind of wild, all Rxs are -1 by defaut?

#okay, let's try re-parameteretizing and just seeing where we get with r.ids
rx_conf = dict(zip(S94_s90_GPR_file.abb_Rx_name, S94_s90_GPR_file.solved_Rx_rule)) #use the original Rx-based ruleset

opt_test2 = CORDA(recon3, rx_conf)
#I was hoping to get out of building the model, but apparently that's not an option...
opt_test2.build()

print(opt_test2)
print([opt_test2.model.reactions.get_by_id(k).reaction for k, used in opt_test2.included.items() if used])
#this prints the list of reactions in the reconstructed model without specifically creating the new model
len(opt_test2.model.reactions)

opt_test2_model = opt_test2.model
for r in opt_test2_model.reactions:
    print(rx_conf[r.id])
#OKAY GOOD, now we have confidences assigned properly. Note that we do have to specify rx_conf, bc we called conf "rx_conf' in this command: opt_test2 = CORDA(recon3, rx_conf)

#which Rxs are high confidence?
for r in opt_test2_model.reactions:
    if rx_conf[r.id] == 3:
        print(r.id)
        print(r)
#biomass is not considered high confidence (i GUESS?)- so what is it?
for r in opt_test2_model.reactions:
    if r.id== "biomass_maintenance":
        print(rx_conf[r.id])
#WTF... my biomass Rx confidence is set to 0!?!?!?!??! No wonder the model doesn't run.
#all right, time to fix THAT. 
rx_conf["biomass_maintenance"] = 3
opt_test3 = CORDA(recon3, rx_conf)
opt_test3_model = opt_test3.model

for r in opt_test3_model.reactions:
    if r.id== "biomass_maintenance":
        print(rx_conf[r.id])
#good stuff. 

opt_test3.build()
#built!

for r in opt_test2_model.reactions:
    if r.id== "biomass_maintenance":
        print(rx_conf[r.id])

#OKAY, PSA: the .build() option seems to be editing recon3 itself, and all Rxs associated with it. I'm... not sure how to handle that. 

#and now, we run the model.
print(opt_test3_model.summary())
#nope, didn't work.let's try setting biomass as objective.
opt_test3_model.objective="biomass_maintenance"
print(opt_test3_model.summary())

#still nothing. The only other thing I can think of to do at this point is to edit all the non-gene-associated Rxs to 1 or 2.
#the rationale here is that, if Rxs were included in the model in the first place, they must be necessary for the model to run. Right?

#done. Let's give it a shot.
#import a test GPR file
S94_s90_GPR_file_V2= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/s94_S90_GPRs_binned_solved_V2.csv")
S94_s90_GPR_file_V2.drop(columns={"Unnamed: 0", "edited_Rx_rule"}, inplace=True)
#also (maybe?) need to split Rx_name column
abb_Rx_name = []
for row in S94_s90_GPR_file_V2["Rx_name"]:
    splitlist = row.split(r": ") #remove everything after :
    newrow = splitlist[0]
    abb_Rx_name.append(newrow)
del row, splitlist, newrow
S94_s90_GPR_file_V2["abb_Rx_name"] = abb_Rx_name

S94_s90_GPR_file_V2.drop(columns=['Rx_name'], inplace=True)     

##########ok, convert GPR file to dictionary
gene_conf_V2 = dict(zip(S94_s90_GPR_file_V2.abb_Rx_name, S94_s90_GPR_file_V2.solved_Rx_rule))

opt_test4 = CORDA(recon3, gene_conf_V2)
opt_test4_model = opt_test4.model

for r in opt_test4_model.reactions:
    if r.id== "biomass_maintenance":
        print(gene_conf_V2[r.id])
#okay, biomass is now confidence=2, which is a good sign. We still want it to be 3, but that's okay. 

gene_conf_V2["biomass_maintenance"] = 3
opt_test4 = CORDA(recon3, gene_conf_V2)
opt_test4_model = opt_test4.model
opt_test4.build()
#built!

opt_test4_model.objective="biomass_maintenance"

print(opt_test4_model.summary())

##########################AAAAAAAA FUCK IT STARTINGOVER
import os
import corda
import pandas as pd
import copy

#first, switch to our NEW FOLDER
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models")

from corda import reaction_confidence
from cobra.io import load_matlab_model
from cobra.io import read_sbml_model
from corda import CORDA    


#import our corda model
recon3 = load_matlab_model("/Users/adamo010/Documents/CORDA_CRC_human_models/Recon3D_301/Recon3DModel_301.mat")
#recon3 is the generic human model that we'll be parameterizing.

#first, let's make sure recon3 isn't irrevocably fucked.
print(recon3.summary())
#it is okay. Copy for safety's sake. 
recon3_edited = copy.deepcopy(recon3)

#import and fix ourfile so it makes an Rx/rule dictionary called gene_conf_V2
S94_s90_GPR_file_V2= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/s94_S90_GPRs_binned_solved_V2.csv")
S94_s90_GPR_file_V2.drop(columns={"Unnamed: 0", "edited_Rx_rule"}, inplace=True)
abb_Rx_name = []
for row in S94_s90_GPR_file_V2["Rx_name"]:
    splitlist = row.split(r": ") #remove everything after :
    newrow = splitlist[0]
    abb_Rx_name.append(newrow)
del row, splitlist, newrow
S94_s90_GPR_file_V2["abb_Rx_name"] = abb_Rx_name
S94_s90_GPR_file_V2.drop(columns=['Rx_name'], inplace=True)     
gene_conf_V2 = dict(zip(S94_s90_GPR_file_V2.abb_Rx_name, S94_s90_GPR_file_V2.solved_Rx_rule))

#set biomass_maintenance as bin=3
gene_conf_V2["biomass_maintenance"] = 3

#CORDA-ize the model
opt_test4A = CORDA(recon3_edited, gene_conf_V2)
opt_test4A.build() #2:58pm

#AND NOW WE WAIT
#checked back in; done at 5:33pm

opt_test4A_model = opt_test4A.model
opt_test4A_model.objective="biomass_maintenance"
print(opt_test4A_model.summary())

#whyyyyyyyyyyyyy
print(opt_test4A)

#I won't know if it's bad that we've only included 181 of 10600 reactions? Why have so many been excluded?
#hmmmmm. Let's change all our -1 reactions to 0s in 01.20.22_combining_GPRs_and_expression_bins_V2.py
#let's try it.
S94_s90_GPR_file_V3= pd.read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_expression_bins_added/s94_S90_GPRs_binned_solved_V2.csv")
S94_s90_GPR_file_V3.drop(columns={"Unnamed: 0", "edited_Rx_rule"}, inplace=True)
abb_Rx_name = []
for row in S94_s90_GPR_file_V3["Rx_name"]:
    splitlist = row.split(r": ") #remove everything after :
    newrow = splitlist[0]
    abb_Rx_name.append(newrow)
del row, splitlist, newrow
S94_s90_GPR_file_V3["abb_Rx_name"] = abb_Rx_name
S94_s90_GPR_file_V3.drop(columns=['Rx_name'], inplace=True)     
gene_conf_V3 = dict(zip(S94_s90_GPR_file_V3.abb_Rx_name, S94_s90_GPR_file_V3.solved_Rx_rule))

#set biomass_maintenance as bin=3
gene_conf_V3["biomass_maintenance"] = 3

#CORDA-ize the model
recon3_edited2 = copy.deepcopy(recon3)
opt_test5 = CORDA(recon3_edited2, gene_conf_V3)
opt_test5.build() #bgin 8:16pm. Finished around 9:50ish the next day

#AND NOW WE WAIT

opt_test5_model = opt_test5.model
opt_test5_model.objective="biomass_maintenance"
print(opt_test5_model.summary())
#STILL ZILCHO
print(opt_test5)
#so, that's better: 271 outof 10600 reactions included. What the hell is happening to the thousands of low/medium reactions?

#going to give this another shot. Now, all non-gene-associated reactions have a bin of 3. 
