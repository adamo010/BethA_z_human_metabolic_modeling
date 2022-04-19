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

#generated these data in 02.11.22_Burns2015_human_models_highfiber_combining_fluxes.py

#for this, borrow from both 12.20.21_Burns2015_Fuso_fluxpath_analysis_exchanges_only_newfilter2.R, and 02.14.22_Burns2022_CORDA_model_analyzing_highfiber_fluxes.R

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
rm(fluxes_unique) #cleanup

##########Step 3: remove all fluxes that aren't exchanges. 
analyzed_fluxes_exchanges <- dplyr::filter(analyzed_fluxes, grepl("EX_",fluxpath_name))
#interesting- no metadata for the exchange reactions. I'll probably have to import the metabolites dataframe
metabolites_metadata <- read.csv("/Users/adamo010/Documents/CORDA_CRC_human_models/recon_model_metabolites.tsv", sep = '\t', header = TRUE)
metabolites_metadata_slimmed <- subset(metabolites_metadata, select= c("abbreviation", "fullName"))
metabolites_metadata_slimmed <-  metabolites_metadata_slimmed %>% rename(fluxpath_name= abbreviation,
                                                            fluxpath_description = fullName)
rm(metabolites_metadata)
#now, remove EX_ substring.
analyzed_fluxes_exchanges <- analyzed_fluxes_exchanges %>%
  mutate_at("fluxpath_name", str_replace, "EX_", "")
#drop blank columns
analyzed_fluxes_exchanges<- subset(analyzed_fluxes_exchanges, select= -c(fluxpath_description, fluxpath_subsystem))
analyzed_fluxes_exchanges_plus_metadata <-left_join(analyzed_fluxes_exchanges, metabolites_metadata_slimmed,
                                                    by = c("fluxpath_name" = "fluxpath_name"))
#how many unique fluxes are we including here? 
ex_fluxes_unique <- unique(analyzed_fluxes_exchanges_plus_metadata$fluxpath_name) #1559 fluxes
rm(ex_fluxes_unique) #cleanup

#########Step 4: filter out zeroes
analyzed_fluxpaths_wide <- select(analyzed_fluxes_exchanges_plus_metadata, -c(SampleID)) #drop unnecessary columns
#first, convert to wide form
analyzed_fluxpaths_wide <- analyzed_fluxpaths_wide %>% 
  spread(Description, flux_value)
#then, drop all rows where both tumor and normal values are nan.
analyzed_fluxpaths_wide_noNANs <- analyzed_fluxpaths_wide[!(is.na(analyzed_fluxpaths_wide$tumor) & is.na(analyzed_fluxpaths_wide$normal)),]
fluxpaths_unique_zerofree <- unique(analyzed_fluxpaths_wide_noNANs$fluxpath_name) #count the number of fluxpaths: still 1559
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
fluxpaths_zerofree_unique <- unique(analyzed_fluxpaths_wide_noNANs_zerofree$fluxpath_name) #165 fluxpaths
#clean up
rm(fluxpaths_unique_zerofree, fluxpaths_zerofree_unique, analyzed_fluxpaths_wide_noNANs, analyzed_fluxpaths_wide)

#########Step 5: filter by flux subsystems which are active in at least half the total number of samples 
#first, create a filtering dataframe that contains counts of each fluxpath_name incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
analyzed_fluxpaths_filtering_df <- analyzed_fluxpaths_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[4:5]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 3 and 4 (tumor and normal) for each row
  group_by(fluxpath_name) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
analyzed_fluxpaths_filter_added <- dplyr::inner_join(analyzed_fluxpaths_wide_noNANs_zerofree, analyzed_fluxpaths_filtering_df, by= "fluxpath_name")
#now, remove all rows where n (the name of the count column) is less than ...44, which is half the samples
analyzed_fluxpaths_filtered <- filter(analyzed_fluxpaths_filter_added, n >= 44)
#count number of unique metabolites left
fluxpaths_filtered_unique <- unique(analyzed_fluxpaths_filtered$fluxpath_name)
#trimming to active in 44 samples gets us to 36 subsystems; 88 gets us to 2 exchanges. Let's try that for now. 
#clean up
rm(fluxpaths_filtered_unique, analyzed_fluxpaths_wide_noNANs_zerofree)

##########Step 6: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe
#drop a couple of columns
analyzed_fluxpaths_paired_only <- select(analyzed_fluxpaths_filtered, -c(difference, n, fluxpath_description))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
#analyzed_metabolites_paired_only <- analyzed_metabolites_paired_only %>% rename(normal= mean_genus_GR_normal, tumor= mean_genus_GR_tumor) #rename these columns
#convert to long form
analyzed_fluxpaths_paired_long <- reshape2::melt(data= analyzed_fluxpaths_paired_only,
                                                 id.vars= c("Patient_Blind_ID", "fluxpath_name"),
                                                 variable.name = "Description",
                                                 value.name = "fluxpath_activity")
#split up by different subsystems and see where we get. 
split_analyzed_fluxpaths_by_flux <- split(analyzed_fluxpaths_paired_long, with(analyzed_fluxpaths_paired_long, interaction(fluxpath_name)), drop = TRUE)
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
  new_value = elem$fluxpath_name[1]
  flux_ids <- c(flux_ids, new_value)}
#merge all lists together. 
fluxpath_statistics <- data.frame(flux_ids, by_flux_pvals, by_flux_qvals)

#########Step 8: merge descriptions back in
metabolites_metadata_slimmed
fluxpath_statistics_plus_metadata <-left_join(fluxpath_statistics, metabolites_metadata_slimmed,
                                                    by = c("flux_ids" = "fluxpath_name"))
##########Step 9: save results as a csv file
fluxpath_statistics_plus_metadata <- fluxpath_statistics_plus_metadata[order(fluxpath_statistics_plus_metadata$by_flux_qvals),] #re-order by q-value
fluxpath_statistics_plus_metadata <- fluxpath_statistics_plus_metadata %>% rename("fluxpath_name"="flux_ids")
write.csv(as.data.frame(fluxpath_statistics_plus_metadata), file="02.17.22_CORDA_flux_exchanges_stats_Burns2015_data_EUstd.csv")
write.csv(as.data.frame(analyzed_fluxpaths_paired_only), file="02.17.22_CORDA_flux_exchanges_alldata_Burns2015_data_EUstd.csv")

##########Step 10: prepare data for graphing START HERE
#LAAAAAAAAAAAAME now I have to calculate differences.
#use analyzed_metabolites_paired_only. First, drop boring cols
#add a new column called diff. tumor-normal
analyzed_fluxpaths_paired_wide = analyzed_fluxpaths_paired_only
analyzed_fluxpaths_paired_wide$diff <- (analyzed_fluxpaths_paired_wide$tumor - analyzed_fluxpaths_paired_wide$normal)
analyzed_fluxpaths_paired_wide <-left_join(analyzed_fluxpaths_paired_wide, metabolites_metadata_slimmed,
                                              by = c("fluxpath_name" = "fluxpath_name"))

#NOW we can add a new column with useful binning of differences for graphing purposes
analyzed_fluxpaths_paired_wide <- analyzed_fluxpaths_paired_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",
  ))

#step 10: NOW GRAPH??!!!
#D-glucose exchange
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("glc_D")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("D-glucose exchange flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.00153, q= 0.0413")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0,1000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Uridine diphosphate glucuronic acid
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("udpglcur")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("UDP glucuronic acid exchange flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.00230, q= 0.0413")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-1000,1000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Nicotinic acid mononucleotide
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("nicrnt")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Nicotinic acid mononucleotide exchange flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.00857, q= 0.103")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-500,1200)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#	Adenosine triphosphate
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("atp")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("ATP exchange flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.0218, q= 0.174")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-500,1200)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Uridine
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("uri")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Uridine exchange flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.0260, q= 0.174")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(0,1200)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Phosphatidylglycerol
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("pglyc_hs")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Phosphatidylglycerol exchange flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.0294, q= 0.174")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-5,0)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Dextrin
ggpaired(subset(analyzed_fluxpaths_paired_wide, fluxpath_name %in% c("dxtrn")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("Dextrin exchange flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.0339, q= 0.174")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average flux,\nmmol/(gDW * h)", limits=c(-100,500)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
