library("ggplot2")
library("readxl")
library("dplyr")
library("tidyverse")
library("ggpubr")
library("ape")
library("lme4")
library("gplots")
library("plotly")
library("tidyr")
library("vegan")
library("data.table")
library("stringr")

#the goal is to graph uptakes of sample-specific CORDA models to see if we can see CRC vs healthy site differences
#there are three datasets to look at; uptake, secretion, and fluxes. This file is for uptake.

#generated these data in 02.11.22_Burns2015_human_models_highfiber_combining_uptakes.py

#starting point is 01.10.22_Niccolai2020_Fuso_fluxpath_analysis_at_subsystem_level.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/CORDA_CRC_human_models/highfiber_diet_model_outputs/") 

##########Step 2: import data and adjust factors
analyzed_uptakes <- read_csv("02.14.22_Burns2015_human_models_highfiber_uptakes_combined.csv") 
analyzed_uptakes$Patient_Blind_ID <- as.factor(analyzed_uptakes$Patient_Blind_ID)
#clean up extra columns
analyzed_uptakes <- rename(analyzed_uptakes, c("C_number"= "C-number", "C_flux"= "C-Flux")) #have to rename to be able to delete; R doesn't like - in column names
analyzed_uptakes <- subset(analyzed_uptakes, select = -c(...1, C_number, abbreviation, Tissue_Tube_ID, Age, Sex, Site, MSI_status, Stage, sample_id))

#how many unique fluxes are we including here? 
uptakes_unique <- unique(analyzed_uptakes$Reaction) #85 uptake reactions
rm(uptakes_unique) #cleanup

########Step 3: filter to only include paired uptakes (i.e uptakes which appear in both tumor and normal samples)
#first, do some reordering- want Patient_Blind_id is next to reaction
#DON't need to resort with this dataset- reaction is already next to Patient_Blind_ID
#analyzed_uptakes <- analyzed_uptakes[, c(1,6,2,3,4,5,7,8)] #want Reaction next to Patient_Blind_ID (NOT metabolite, as there are often sinks/exports of same metab)
#then, create a new dataframe where only metabolites with activity in BOTH tumor and normal samples (i.e. Patient_Blind_id appears twice) are included
analyzed_uptakes_paired <- analyzed_uptakes %>% unite("sorting_col", Reaction:Patient_Blind_ID, remove = FALSE)
analyzed_uptakes_paired_only <- subset(analyzed_uptakes_paired,duplicated(sorting_col) | duplicated(sorting_col, fromLast=TRUE))

#now how many unique fluxes do we get?
paired_uptakes_unique <- unique(analyzed_uptakes_paired_only$Reaction) #26 uptake reactions
rm(paired_uptakes_unique) #cleanup

#MIGHT need to filter by metabs which appear in at least half the samples; stay tuned

#########Step 4: first pass: split into multiple dataframes by metabolite and do stats on each metabolite
split_analyzed_uptakes_by_metab <- split(analyzed_uptakes_paired_only, with(analyzed_uptakes_paired_only, interaction(Reaction)), drop = TRUE)
by_metab_stats <- lapply(split_analyzed_uptakes_by_metab, function(df){wilcox.test(Flux~Description, data=df, exact= FALSE, paired= TRUE)})

##########Step 5: statistics
#make lists of p values
by_metab_pvals <- c() #initialize list
for (elem in by_metab_stats){
  new_value = elem$p.value
  by_metab_pvals <- c(by_metab_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_metab_qvals <- p.adjust(by_metab_pvals, method = "fdr")
#make a list of genus_ids
reaction_ids <- c()
for(elem in split_analyzed_uptakes_by_metab){
  new_value = elem$Reaction[1]
  reaction_ids <- c(reaction_ids, new_value)}
#merge all lists together. 
metabolite_statistics <- data.frame(reaction_ids, by_metab_pvals, by_metab_qvals)

#ALL 1s or NAN. When present, all metabolites have the same value across all samples. Cool



