#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 14:43:09 2021

@author: adamo010
"""
#all right, we're finally going to get to it; metabolic modeling of human host tissue. Thanks to Osbaldo Resendis for writing a python wrapper
#for a matlab version of CORDA: Cost Optimization Reaction Dependency Assessment. This is from Amina Qutub's lab, and is many of the options
#for creating tissue-specific metabolic models. I am using it because a) the paper made something that is most similar to the kinds of models
#I want to generate, and b) b/c I saw Dr. Qutub speak at a conference and she seemed nice.

#we're working in a new virutal environment, CORDA_human_fba, for this, so it might be slow.
import os
import corda

#first, switch to our NEW FOLDER
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models")

from corda import reaction_confidence

gene_conf = {"gene1": 1, "gene2": 3, "gene4": -1} # missing entries are automatically assigned zeroes
rule = "gene1 and gene2 or (gene3 and gene4)"

reaction_confidence(rule, gene_conf)
#thrilling. 

#performing a small reconstruction
from corda import test_model #this is for central carbon metabolism

mod = test_model()
len(mod.reactions) #gives the number of reactions in the model; it is 60.

mod.reactions[59].reaction #the last reaction is the biomass reaction. prints out a nice Rx

#We can now reconstruct the smallest model that still allows for the production of biomass. 
#Let's start by setting the biomass reaction as high confidence and all other reactions as "not include".

conf = {}
for r in mod.reactions: conf[r.id] = -1 #I assume this means only include the last reaction
conf["r60"] = 3

#okay, so conf looks like it's a dictionary of all reactions in the model. each has a value of -1 (from line 40)
#except the 60th element ("r60", which is the biomass reaction), which has a confidence value of 3.

#now, reconstruct the model:

from corda import CORDA

opt = CORDA(mod, conf) #create a new model, opt, from mod(the original model) with confidence values from conf
opt.build() #build the model
print(opt)
    
#speedy as all fuck, I'll say that much.     

#all right, take a look at the actual model:
    
print([opt.model.reactions.get_by_id(k).reaction for k, used in opt.included.items() if used])    
#prints a nice list of all the reactions used in opt.

#We can also define additional metabolic functions that should be possible in the reconstruction. 
#Let's assume we want to be able to produce pep.

opt = CORDA(mod, conf, met_prod="pep")
opt.build()
print(opt)

#presumably pep is in the model already?
#my guess is yes. Pep is a metabolite name present in the model. 

rec = opt.cobra_model("plus_pep")
use = rec.metabolites.pep.reactions
print("# of redundant pathway for pep =", opt.redundancies["EX_CORDA_0"])
for r in use: print(r.reaction)

#there are two Rxs that produce pep. By default CORDA uses redundancy. This means, in case there are several 
#minimal pathways to reach your objective, CORDA will include all of those (which is good since it gives your 
#model some robustness). In this case this did happen and yields two different forms to produce pep. 
#If we do not want that feature we can modify the parameter n in the CORDA initializer which denotes the 
#maximum number of redundant pathways to include.

#to only include a single pathway, do the following:
opt = CORDA(mod, conf, met_prod="pep", n=1)
opt.build()

rec_min = opt.cobra_model("plus_pep_nored")
print("used", len(rec_min.reactions), "reactions")
print("# of redundant pathway for pep =", opt.redundancies["EX_CORDA_0"])
use = rec_min.metabolites.pep.reactions
for r in use: print(r.reaction)

#nice. Now, to figure out how to scale this to 20,000 host genes...