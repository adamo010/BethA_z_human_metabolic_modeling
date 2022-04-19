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
os.chdir("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_analyzed_data")

expression_data = pd.read_csv('all_samples_protein_coding_subread_counts.txt', delimiter = "\t")
expression_data.rename(columns={"Unnamed: 0": "gene_id"}, inplace=True)

#now, plot data: can I make 10 plots of randomly selected sets of 100 rows?
expression_data.plot.bar()



