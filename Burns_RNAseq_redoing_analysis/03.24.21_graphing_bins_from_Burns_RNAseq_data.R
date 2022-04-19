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
library("ggpubr")
library("rstatix")

setwd("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/TPM_abundances_from_kallisto")

counts_df<- read_csv("03.24.21_counts_of_binned_RNAseq_data.csv", col_names = TRUE)

#set first column (bins) as factor
#counts_df$X1 <- as.factor(counts_df$X1)

#transpose dataframe and save as new dataframe
counts_df_flip <- as.data.frame(t(counts_df))

#convert first row to column names
colnames(counts_df_flip) <- as.character(unlist(counts_df_flip[1,]))
counts_df_flip = counts_df_flip[-1, ]

#I would also like to add a Description column, to see if there are differences between CRC and normal cell expression patterns

#merge dataframes
#Assss the column names in counts_df_flip don't match that in the metadata.
#copy the column containing row_names, trim off everything after the underscore, and rename it as "tube_id"
counts_df_flip2 <- counts_df_flip
counts_df_flip2 <- tibble::rownames_to_column(counts_df_flip2, "tube_id") 
counts_df_flip2$tube_id <- sub("_[^_]+$", "", counts_df_flip2$tube_id)

#import metadata
metadata <- read.table(file = "/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_analyzed_data/CRC_metadata_rnaseq.txt", sep = "\t", header=TRUE)

#joined_tibble <- left_join(mytibble, mylookup_tibble, by = c("OP_UNIQUE_CARRIER" = "Code"))

counts_with_metadata <- left_join(counts_df_flip2, metadata, by= c("tube_id" = "Tissue.RNA.DNA_Tube_ID"))

#plot: I guess a scatterplot
#UGH I'm going to have to convert this to long form, aren't I.
colnames(counts_with_metadata)[2] <- "0_(not_inlcuded)"
colnames(counts_with_metadata)[3] <- "1_(low_expression)"
colnames(counts_with_metadata)[4] <- "-1_(no_expression)"
colnames(counts_with_metadata)[5] <- "2_(medium_expression)"
colnames(counts_with_metadata)[6] <- "3_(high_expression)"

counts_with_metadata_long <- gather(counts_with_metadata, expression_bin, gene_counts,
                                    "0_(not_inlcuded)":"3_(high_expression)", factor_key=TRUE)

#okay NOW do a scatterplot.
#y values are gene_counts
#x values are expression_bin
#each dot is a Patient_Blind_ID
#facet_wrap by Description

ggplot(counts_with_metadata_long, aes(expression_bin, gene_counts)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               stackdir ="center",
               dotsize = 0.5,
               fill= "blue") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6))+
  facet_wrap(~ Description)

#not a bad first try

#let's futz with it now

#first, change the x axis order to be sensible
counts_with_metadata_long$expression_bin <- factor(counts_with_metadata_long$expression_bin,
                                                   levels = c("0_(not_inlcuded)","-1_(no_expression)","1_(low_expression)",
                                                              "2_(medium_expression)","3_(high_expression)"))

#first, plot each patient separately to compare tumor vs normal counts
ggplot(counts_with_metadata_long, aes(expression_bin, gene_counts, fill=factor(Description), color=factor(Description))) +
  #geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               stackdir ="center",
               dotsize = 0.5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6), plot.title=element_text(hjust=0.5))+
  labs(title="Number of genes in each expression bin,\nby Patient_Blind_ID")+ 
  facet_wrap(~ Patient_Blind_ID)

#this is actually not that helpful, let's color by Patient_Blind_ID and facet wrap by description
ggplot(counts_with_metadata_long, aes(expression_bin, gene_counts, fill=factor(Patient_Blind_ID), color=factor(Patient_Blind_ID))) +
  #geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               stackdir ="center",
               dotsize = 0.5) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6), plot.title=element_text(hjust=0.5))+
  labs(title="Number of genes in each expression bin,\nby Patient_Blind_ID")+ 
  facet_wrap(~ Description)

#also not that interesting. 
#don't bother with dotplot, just graph boxplot
ggplot(counts_with_metadata_long, aes(expression_bin, gene_counts)) +
  geom_boxplot() +
  #geom_dotplot(binaxis = 'y',
              # stackdir ="center",
               #dotsize = 0.5) +
  theme(axis.text.x = element_text(angle=40, vjust = 1, hjust=0.9), plot.title=element_text(hjust=0.5))+
  labs(title="Number of genes in each expression bin,\nby Patient_Blind_ID")+ 
  facet_wrap(~ Description)


#let's back up. what do I actually care about?
#do the counts of genes in each category differ between CRC and normal samples (paired?)
#in which case, will want to facet_wrap by expression bin

ggplot(counts_with_metadata_long, aes(Description, gene_counts, fill=factor(Description), color=factor(Description))) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir ="center", dotsize = 0.5) +
  theme(axis.text.x = element_text(angle=40, vjust = 1, hjust=0.9), plot.title=element_text(hjust=0.5), legend.position="none")+
  labs(title="Number of genes in each expression bin,\nby tumor/normal status")+ 
  facet_wrap(~ expression_bin, scales= "free")

#just do boxplots

ggplot(counts_with_metadata_long, aes(Description, gene_counts, color=factor(Description))) +
  geom_dotplot(binaxis = 'y', stackdir ="center", dotsize = 0.5) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=40, vjust = 1, hjust=0.9), plot.title=element_text(hjust=0.5), legend.position="none")+
  labs(title="Number of genes in each expression bin,\nby tumor/normal status") + 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  facet_wrap(~ expression_bin, scales= "free")

#these stats are inscrutable. levels = c("0_(not_inlcuded)","-1_(no_expression)","1_(low_expression)",
#"2_(medium_expression)","3_(high_expression)"))

#another way of doing this:(and the correct one)
########################Level 0: Not included########################
#subset data and run stats on it
bin0stat <- subset(counts_with_metadata_long, expression_bin =="0_(not_inlcuded)", select= c(Description, gene_counts))
wilcox.test(gene_counts ~ Description, data = bin0stat, correct= FALSE, paired = TRUE)
#p-value = 0.2483

ggplot(subset(counts_with_metadata_long, expression_bin %in% c("0_(not_inlcuded)")), aes(Description, gene_counts, color= Description)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir ="center", dotsize = 0.5, aes(fill= Description)) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  labs(title="Number of genes in bin 0 (absent),\nby tumor/normal status, p=0.2483") +
  theme_linedraw() + theme_light() +
  theme(plot.title=element_text(hjust=0.5), legend.position="none") 

########################Level -1: No expression########################
#subset data and run stats on it
binneg1stat <- subset(counts_with_metadata_long, expression_bin =="-1_(no_expression)", select= c(Description, gene_counts))
wilcox.test(gene_counts ~ Description, data = binneg1stat, correct= FALSE, paired = TRUE)
#p-value = 0.8702

ggplot(subset(counts_with_metadata_long, expression_bin %in% c("-1_(no_expression)")), aes(Description, gene_counts, color= Description)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir ="center", dotsize = 0.5, aes(fill= Description)) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  labs(title="Number of genes in bin -1 (no expression),\nby tumor/normal status, p=0.8702") +
  theme_linedraw() + theme_light() +
  theme(plot.title=element_text(hjust=0.5), legend.position="none") 

########################Level 1: Low expression########################
#subset data and run stats on it
bin1_stat <- subset(counts_with_metadata_long, expression_bin =="1_(low_expression)", select= c(Description, gene_counts))
wilcox.test(gene_counts ~ Description, data = bin1_stat, correct= FALSE, paired = TRUE)
#p-value = 0.8308

ggplot(subset(counts_with_metadata_long, expression_bin %in% c("1_(low_expression)")), aes(Description, gene_counts, color= Description)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir ="center", dotsize = 0.5, aes(fill= Description)) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  labs(title="Number of genes in bin 1 (low expression),\nby tumor/normal status, p=0.8308") +
  theme_linedraw() + theme_light() +
  theme(plot.title=element_text(hjust=0.5), legend.position="none") 

########################Level 2: Medium expression########################
#subset data and run stats on it
bin2_stat <- subset(counts_with_metadata_long, expression_bin =="2_(medium_expression)", select= c(Description, gene_counts))
wilcox.test(gene_counts ~ Description, data = bin2_stat, correct= FALSE, paired = TRUE)
#p-value = 0.007586

ggplot(subset(counts_with_metadata_long, expression_bin %in% c("2_(medium_expression)")), aes(Description, gene_counts, color= Description)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir ="center", dotsize = 0.5, aes(fill= Description)) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  labs(title="Number of genes in bin 2 (medium expression),\nby tumor/normal status, p=0.007586") +
  theme_linedraw() + theme_light() +
  theme(plot.title=element_text(hjust=0.5), legend.position="none") 

########################Level 3: High expression########################
#subset data and run stats on it
bin3_stat <- subset(counts_with_metadata_long, expression_bin =="3_(high_expression)", select= c(Description, gene_counts))
wilcox.test(gene_counts ~ Description, data = bin3_stat, correct= FALSE, paired = TRUE)
#p-value = 0.008758
ggplot(subset(counts_with_metadata_long, expression_bin %in% c("3_(high_expression)")), aes(Description, gene_counts, color= Description)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir ="center", dotsize = 0.5, aes(fill= Description)) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  labs(title="Number of genes in bin 3 (high expression),\nby tumor/normal status, p=0.008758") +
  theme_linedraw() + theme_light() +
  theme(plot.title=element_text(hjust=0.5), legend.position="none") 



