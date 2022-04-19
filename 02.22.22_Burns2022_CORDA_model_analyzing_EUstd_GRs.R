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
#there are three datasets to look at; uptake, secretion, and fluxes. This file is for fluxes.

#generated these data in 02.16.22_Burns2015_human_models_EUstd_combining_fluxes.py

#starting point is 02.14.22_Burns2022_CORDA_model_analyzing_highfiber_fluxes_exchanges_only.R

##########Preliminary steps
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/CORDA_CRC_human_models/EUstd_diet_model_outputs/") 

##########Step 2: import data and adjust factors
group_GR_all <- read.csv(file= "02.22.22_Burns2015_human_models_EUstd_GRs_combined.csv", header=FALSE)
group_GR_all <- group_GR_all %>% rename(sample_id = V1,
                                        growth_rate = V2)
#import metadata
metadata= read_csv("/Users/adamo010/Documents/CORDA_CRC_human_models/CRC_metadata_rnaseq_edited_by_BA.csv")
metadata_slimmed <- select(metadata, c(Tissue_Tube_ID, Patient_Blind_ID, SampleID, Description)) #drop the columns called ...1 (the index column, and the sample ID column)

#merge metadata with data
analyzed_GRs_plus_metadata <-left_join(group_GR_all, metadata_slimmed, by = c("sample_id" = "Tissue_Tube_ID"))

#set factors
analyzed_GRs_plus_metadata$Patient_Blind_ID <- as.factor(analyzed_GRs_plus_metadata$Patient_Blind_ID)
analyzed_GRs_plus_metadata <- select(analyzed_GRs_plus_metadata, -c(SampleID, sample_id))

###########Step 3: rearranging
#split GR column into 2 columns, named by different values of Description
analysed_GRs_wide = analyzed_GRs_plus_metadata %>% spread(Description, growth_rate) 
#now, create a new column called mean_GR_difference
analysed_GRs_wide <- analysed_GRs_wide %>%
  mutate(GR_difference=tumor-normal)
#now, create another new column called Difference_direction
analysed_GRs_wide <- analysed_GRs_wide %>%
  mutate(Difference_direction = case_when(
    GR_difference > 0 ~ "tumor > normal",
    GR_difference < 0 ~ "tumor < normal",
    GR_difference == 0 ~ "Zero change",
  ))

############Step 4: stats and graphing
#cool. Now let's make some graphs and run some stats.
# Compute t-test- REMEMBER to use the non wide form of the dataframe
#GR_by_description <- wilcox.test(growth_rate ~ Description, data = analyzed_GRs_plus_metadata, paired = TRUE)
#can't use wilcoxan- there are pairs with differences=0, so it won't work. 
GR_by_description <- t.test(growth_rate ~ Description, data = analyzed_GRs_plus_metadata, paired = TRUE)
GR_by_description$p.value #print p-value #aha- 0.0205

############Step 5: graph
ggpaired(analysed_GRs_wide, cond1= "normal", cond2= "tumor",
         y= "growth_rate", fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Colonocyte growth rates, Burns2015 data"), subtitle = expression("paired T-test, p=0.0205")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Colonocyte growth\nrate, mmol/(gDW * h)", limits=c(0,300)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill='none')





