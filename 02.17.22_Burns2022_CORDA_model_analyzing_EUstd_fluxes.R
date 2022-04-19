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
analyzed_fluxes <- read_csv("02.16.22_Burns2015_human_models_EUstd_fluxes_combined.csv") 
analyzed_fluxes$Patient_Blind_ID <- as.factor(analyzed_fluxes$Patient_Blind_ID)
#clean up extra columns START HERE.
analyzed_fluxes <- subset(analyzed_fluxes, select = -c(...1, sample_id, Tissue_Tube_ID, Age, Sex, Site, MSI_status, Stage))
#need to do some renaming also. for renaming, newname = old name
analyzed_fluxes <- analyzed_fluxes %>% rename(fluxpath_subsystem= subsystem,
                                                    flux_value = fluxes,
                                                    fluxpath_description = description)
#how many unique fluxes are we including here? 
fluxes_unique <- unique(analyzed_fluxes$fluxpath_name) #6318 fluxes
#how many subsystems?
subsystems_unique <- unique(analyzed_fluxes$fluxpath_subsystem) #97 subsystems
rm(fluxes_unique, subsystems_unique) #cleanup

#########Step 3: collapse this dataset by subsystem; average across fluxpaths within a subsystem within a sample.
#drop unnecessary columns
analyzed_fluxes_subset <- subset(analyzed_fluxes, select = -c(fluxpath_description))
#create a new sorting column to get unique values for each Sample_ID/Fluxpath_subsystem combo
analyzed_fluxes_subset <- analyzed_fluxes_subset %>% 
  mutate(sorting_col = paste0(SampleID, "_", fluxpath_subsystem)) 
#AAAA I think I have to replace NANs with 0s. Otherwise the averages won't work 
#analyzed_fluxpaths_subset$fluxpath_amount[is.na(analyzed_fluxpaths_subset$fluxpath_amount)] = 0
#hmmmm, is it better to convert to zeroes or drop??? Probably drop. Let's go with that.
analyzed_fluxes_subset_noNANs <- analyzed_fluxes_subset[!is.na(analyzed_fluxes_subset$flux_value),]
#then average across sorting col:
analyzed_fluxes_subset_averaged <- analyzed_fluxes_subset_noNANs %>%        # Specify data frame
  group_by(sorting_col, SampleID, Description, Patient_Blind_ID, fluxpath_subsystem) %>%    # Specify group indicator (columns to keep)
  summarise_at(vars(flux_value),                   # Specify column
               list(av_fluxpath_amount = mean)) %>% # Specify function
  ungroup() #do this so select (dropping columns) can happen later
#nice. End up with a dataframe where each row is a unique subsystem/sample average of fluxes from pathways within that subsystem
#how many unique fluxes are we including here? 
fluxes_averaged_unique <- unique(analyzed_fluxes_subset_averaged$fluxpath_subsystem) #97 subsystems- no loss from original data
#clean up
rm(fluxes_averaged_unique, analyzed_fluxes_subset_noNANs, analyzed_fluxes_subset)

#########Step 4: filter out zeroes
analyzed_fluxpaths_wide <- select(analyzed_fluxes_subset_averaged, -c(sorting_col, SampleID)) #drop unnecessary columns
#first, convert to wide form
analyzed_fluxpaths_wide <- analyzed_fluxpaths_wide %>% 
  spread(Description, av_fluxpath_amount)
#then, drop all rows where both tumor and normal values are nan.
analyzed_fluxpaths_wide_noNANs <- analyzed_fluxpaths_wide[!(is.na(analyzed_fluxpaths_wide$tumor) & is.na(analyzed_fluxpaths_wide$normal)),]
fluxpaths_unique_zerofree <- unique(analyzed_fluxpaths_wide_noNANs$fluxpath_subsystem) #count the number of fluxpaths: still 97
#for this dataset, there are no nans and all data are kept
#replace remaining NA values with zeroes
analyzed_fluxpaths_wide_noNANs$tumor[is.na(analyzed_fluxpaths_wide_noNANs$tumor)] = 0
analyzed_fluxpaths_wide_noNANs$normal[is.na(analyzed_fluxpaths_wide_noNANs$normal)] = 0
#create a difference column
analyzed_fluxpaths_wide_noNANs$difference <- (analyzed_fluxpaths_wide_noNANs$tumor - analyzed_fluxpaths_wide_noNANs$normal)
#remove any rows where differences are 0
analyzed_fluxpaths_wide_noNANs_zerofree <- analyzed_fluxpaths_wide_noNANs[(analyzed_fluxpaths_wide_noNANs$difference !=0),]
#NOTE: for this dataset, a few rows are removed but no fluxpaths
#how many unique fluxes are we including here? 
fluxpaths_zerofree_unique <- unique(analyzed_fluxpaths_wide_noNANs_zerofree$fluxpath_subsystem) #73 fluxpaths
#clean up
rm(fluxpaths_unique_zerofree, fluxpaths_zerofree_unique, analyzed_fluxpaths_wide_noNANs, analyzed_fluxpaths_wide, analyzed_fluxes_subset_averaged)

#########Step 5: filter by flux subsystems which are active in at least half the total number of samples 
#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
analyzed_fluxpaths_filtering_df <- analyzed_fluxpaths_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:4]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 3 and 4 (tumor and normal) for each row
  group_by(fluxpath_subsystem) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
analyzed_fluxpaths_filter_added <- dplyr::inner_join(analyzed_fluxpaths_wide_noNANs_zerofree, analyzed_fluxpaths_filtering_df, by= "fluxpath_subsystem")
#now, remove all rows where n (the name of the count column) is less than ...44, which is half the samples
analyzed_fluxpaths_filtered <- filter(analyzed_fluxpaths_filter_added, n >= 44)
#count number of unique metabolites left
fluxpaths_filtered_unique <- unique(analyzed_fluxpaths_filtered$fluxpath_subsystem)
#trimming to active in 44 samples gets us to 39 subsystems; 88 gets us to 11 subsystems. Let's try that for now. 
#clean up
rm(fluxpaths_filtered_unique, analyzed_fluxpaths_wide_noNANs_zerofree)

##########Step 6: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe
#drop a couple of columns
analyzed_fluxpaths_paired_only <- select(analyzed_fluxpaths_filtered, -c(difference, n))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
#analyzed_metabolites_paired_only <- analyzed_metabolites_paired_only %>% rename(normal= mean_genus_GR_normal, tumor= mean_genus_GR_tumor) #rename these columns
#convert to long form
analyzed_fluxpaths_paired_long <- reshape2::melt(data= analyzed_fluxpaths_paired_only,
                                                 id.vars= c("Patient_Blind_ID", "fluxpath_subsystem"),
                                                 variable.name = "Description",
                                                 value.name = "fluxpath_activity")
#split up by different subsystems and see where we get. 
split_analyzed_fluxpaths_by_flux <- split(analyzed_fluxpaths_paired_long, with(analyzed_fluxpaths_paired_long, interaction(fluxpath_subsystem)), drop = TRUE)
by_flux_stats <- lapply(split_analyzed_fluxpaths_by_flux, function(df){wilcox.test(fluxpath_activity~Description, data=df, exact= FALSE, paired= TRUE)})

##########Step 7: statistics
#make lists of p values
by_flux_pvals <- c() #initialize list
for (elem in by_flux_stats){
  new_value = elem$p.value
  by_flux_pvals <- c(by_flux_pvals, new_value)}
rm(new_value)
#make a list of q-values
by_flux_qvals <- p.adjust(by_flux_pvals, method = "fdr")
#make a list of genus_ids
flux_ids <- c()
for(elem in split_analyzed_fluxpaths_by_flux){
  new_value = elem$fluxpath_subsystem[1]
  flux_ids <- c(flux_ids, new_value)}
#merge all lists together. 
fluxpath_statistics <- data.frame(flux_ids, by_flux_pvals, by_flux_qvals)

##########Step 8: save results as a csv file
fluxpath_statistics <- fluxpath_statistics[order(fluxpath_statistics$by_flux_qvals),] #re-order by q-value
fluxpath_statistics <- fluxpath_statistics %>% rename("fluxpath_name"="flux_ids")
write.csv(as.data.frame(fluxpath_statistics), file="02.17.22_CORDA_flux_subsystems_stats_Burns2015_data_EUstd.csv")

##########Step 9: prepare data for graphing START HERE
#LAAAAAAAAAAAAME now I have to calculate differences.
#use analyzed_metabolites_paired_only. First, drop boring cols
#add a new column called diff. tumor-normal
analyzed_fluxpaths_paired_wide = analyzed_fluxpaths_paired_only
analyzed_fluxpaths_paired_wide$diff <- (analyzed_fluxpaths_paired_wide$tumor - analyzed_fluxpaths_paired_wide$normal)

#NOW we can add a new column with useful binning of differences for graphing purposes
analyzed_fluxpaths_paired_wide <- analyzed_fluxpaths_paired_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",
  ))

#step 10: NOW GRAPH??!!!
#Pyrimidine catabolism
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_subsystem %in% c("Pyrimidine catabolism")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Pyrimidine catabolism flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.00333, q= 0.127")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-150,150)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Bile acid synthesis
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_subsystem %in% c("Bile acid synthesis")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Bile acid synthesis fluxes, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.0101, q= 0.191")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0,80)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Glutamate metabolism
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_subsystem %in% c("Glutamate metabolism")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Glutamate metabolism fluxes, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.0362, q= 0.267")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-200,300)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Lysine metabolism
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_subsystem %in% c("Lysine metabolism")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Lysine metabolism fluxes, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.0247, q= 0.267")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-150,50)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

