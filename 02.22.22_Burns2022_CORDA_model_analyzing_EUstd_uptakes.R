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
analyzed_uptakes <- read_csv("02.16.22_Burns2015_human_models_EUstd_uptakes_combined.csv") 
analyzed_uptakes$Patient_Blind_ID <- as.factor(analyzed_uptakes$Patient_Blind_ID)
#clean up extra columns START HERE.
analyzed_uptakes <- subset(analyzed_uptakes, select = c(Metabolite, Reaction, Flux, sample_id, fullName, Patient_Blind_ID, Description))

#need to do some renaming also. for renaming, newname = old name
analyzed_uptakes <- analyzed_uptakes %>% rename(flux_value = Flux)
#how many unique fluxes are we including here? 
uptakes_unique <- unique(analyzed_uptakes$Reaction) #85 uptake rxs
rm(uptakes_unique) #cleanup

#########Step 3: filter out zeroes
analyzed_uptakes_wide <- select(analyzed_uptakes, -c(sample_id, Metabolite, fullName)) #drop unnecessary columns
#first, convert to wide form
#hmmmm, I have some metabolites that don't appear in paired samples. how do I... handle this?
#have to keep the Reaction column, b/c there are different values in this column vs in metabolites (e.g. exported vs sink proline)
analyzed_uptakes_wide <- analyzed_uptakes_wide %>% 
  spread(Description, flux_value)
#then, drop all rows where both tumor and normal values are nan.
analyzed_uptakes_wide_noNANs <- analyzed_uptakes_wide[!(is.na(analyzed_uptakes_wide$tumor) & is.na(analyzed_uptakes_wide$normal)),]
uptakes_unique_zerofree <- unique(analyzed_uptakes_wide_noNANs$Reaction) #count the number of fluxpaths: 85
#for this dataset, there are no nans and all data are kept
#replace remaining NA values with zeroes
analyzed_uptakes_wide_noNANs$tumor[is.na(analyzed_uptakes_wide_noNANs$tumor)] = 0
analyzed_uptakes_wide_noNANs$normal[is.na(analyzed_uptakes_wide_noNANs$normal)] = 0
#create a difference column
analyzed_uptakes_wide_noNANs$difference <- (analyzed_uptakes_wide_noNANs$tumor - analyzed_uptakes_wide_noNANs$normal)
#remove any rows where differences are 0
analyzed_uptakes_wide_noNANs_zerofree <- analyzed_uptakes_wide_noNANs[(analyzed_uptakes_wide_noNANs$difference !=0),]
#NOTE: for this dataset, a few rows are removed but no fluxpaths
#how many unique fluxes are we including here? 
uptakes_zerofree_unique <- unique(analyzed_uptakes_wide_noNANs_zerofree$Reaction) #83 fluxpaths
#clean up
rm(uptakes_unique_zerofree, uptakes_zerofree_unique, analyzed_uptakes_wide_noNANs, analyzed_uptakes_wide)

#########Step 4: filter by uptakes which are active in at least half the total number of samples 
#first, create a filtering dataframe that contains counts of each Reaction incidence in the dataframe
#need to have n=2 if values are nonzero in both tumor and normal columns; n=1 if values are nonzero in only one column
analyzed_uptakes_filtering_df <- analyzed_uptakes_wide_noNANs_zerofree %>%  #create a new dataframe by filtering the old one
  mutate(count_of_nonzeros = rowSums(.[3:4]!=0)) %>% #create a new column, count_of_nonzeroes, that counts # nonzero values in columns 3 and 4 (tumor and normal) for each row
  group_by(Reaction) %>% #group by the column of interest
  summarize(n = sum(count_of_nonzeros)) %>% #sums the count_of_nonzeros values within each fluxpath_subsystem and saves the value as n
  mutate(n) #add n as a new column in the dataframe
#create a new dataframe that merges the old and new dataframes
analyzed_uptakes_filter_added <- dplyr::inner_join(analyzed_uptakes_wide_noNANs_zerofree, analyzed_uptakes_filtering_df, by= "Reaction")
#now, remove all rows where n (the name of the count column) is less than ...44, which is half the samples
analyzed_uptakes_filtered <- filter(analyzed_uptakes_filter_added, n >= 44)
#count number of unique metabolites left
uptakes_filtered_unique <- unique(analyzed_uptakes_filtered$Reaction)
#trimming to active in 44 samples gets us to 37 subsystems; 88 gets us to 2 exchanges. Let's try that for now. 
#clean up
rm(uptakes_filtered_unique, analyzed_uptakes_wide_noNANs_zerofree)

##########Step 5: split into multiple dataframes by metabolite and do stats on each metabolite
#first, edit the dataframe
#drop a couple of columns
analyzed_uptakes_paired_only <- select(analyzed_uptakes_filtered, -c(difference, n))
#rename a couple of columns to make life easier later. Yeah, yeah, I know, this is the reverse of what I did before.
#analyzed_metabolites_paired_only <- analyzed_metabolites_paired_only %>% rename(normal= mean_genus_GR_normal, tumor= mean_genus_GR_tumor) #rename these columns
#convert to long form
analyzed_uptakes_paired_long <- reshape2::melt(data= analyzed_uptakes_paired_only,
                                                 id.vars= c("Patient_Blind_ID", "Reaction"),
                                                 variable.name = "Description",
                                                 value.name = "fluxpath_activity")
#split up by different subsystems and see where we get. 
split_analyzed_uptakes_by_flux <- split(analyzed_uptakes_paired_long, with(analyzed_uptakes_paired_long, interaction(Reaction)), drop = TRUE)
by_flux_stats <- lapply(split_analyzed_uptakes_by_flux, function(df){wilcox.test(fluxpath_activity~Description, data=df, exact= FALSE, paired= TRUE)})

##########Step 6: statistics
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
for(elem in split_analyzed_uptakes_by_flux){
  new_value = elem$Reaction[1]
  flux_ids <- c(flux_ids, new_value)}
#merge all lists together. 
fluxpath_statistics <- data.frame(flux_ids, by_flux_pvals, by_flux_qvals)

#########Step 7: merge descriptions back in
metabolites_metadata_slimmed <- select(analyzed_uptakes, c(Reaction, fullName))
metabolites_metadata_slimmed <- metabolites_metadata_slimmed[!duplicated(metabolites_metadata_slimmed$Reaction), ]
fluxpath_statistics_plus_metadata <-left_join(fluxpath_statistics, metabolites_metadata_slimmed,
                                              by = c("flux_ids" = "Reaction"))
##########Step 8: save results as a csv file
fluxpath_statistics_plus_metadata <- fluxpath_statistics_plus_metadata[order(fluxpath_statistics_plus_metadata$by_flux_qvals),] #re-order by q-value
fluxpath_statistics_plus_metadata <- fluxpath_statistics_plus_metadata %>% rename("Reaction"="flux_ids")
write.csv(as.data.frame(fluxpath_statistics_plus_metadata), file="02.22.22_CORDA_uptake_stats_Burns2015_data_EUstd.csv")

##########Step 9: prepare data for graphing START HERE
#LAAAAAAAAAAAAME now I have to calculate differences.
#use analyzed_metabolites_paired_only. First, drop boring cols
#add a new column called diff. tumor-normal
analyzed_uptakes_paired_wide = analyzed_uptakes_paired_only
analyzed_uptakes_paired_wide$diff <- (analyzed_uptakes_paired_wide$tumor - analyzed_uptakes_paired_wide$normal)
analyzed_uptakes_paired_wide <-left_join(analyzed_uptakes_paired_wide, metabolites_metadata_slimmed,
                                           by = c("Reaction" = "Reaction"))

#NOW we can add a new column with useful binning of differences for graphing purposes
analyzed_uptakes_paired_wide <- analyzed_uptakes_paired_wide %>%
  mutate(Difference_direction = case_when(
    diff > 0 ~ "tumor > normal",
    diff < 0 ~ "tumor < normal",
    diff == 0 ~ "Zero change",
  ))

#step 10: NOW GRAPH??!!!
#L-Proline sink
ggpaired(subset(analyzed_uptakes_paired_wide, Reaction %in% c("sink_pro_L[c]")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("L-Proline sink uptake flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.00144, q= 0.0397")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average uptake flux,\nmmol/(gDW * h)", limits=c(0,1000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#Uridine diphosphate glucuronic acid export
ggpaired(subset(analyzed_uptakes_paired_wide, Reaction %in% c("EX_udpglcur[e]")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("UDP glucuronic acid uptake flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.00214, q= 0.0397")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average uptake flux,\nmmol/(gDW * h)", limits=c(0,1000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")

#L-Glutamine sink
ggpaired(subset(analyzed_uptakes_paired_wide, Reaction %in% c("sink_gln_L[c]")), cond1= "normal", cond2= "tumor",
         fill = "condition", line.color = "Difference_direction") +
  labs(title = expression("L-glutamine sink flux, Burns2015 data (human)"), subtitle = expression("paired Wilcoxan, p=0.00821, q= 0.0799")) +
  theme(plot.title = element_text(hjust=0.5), plot.subtitle= element_text(hjust=0.5)) +
  scale_y_continuous(name = "Average uptake flux,\nmmol/(gDW * h)", limits=c(0,1000)) +
  scale_x_discrete(labels = c('Normal','Tumor')) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.text=element_text(size=11)) +
  guides(fill="none")
