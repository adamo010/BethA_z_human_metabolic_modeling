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

#so, in order to use CORDA to make human metabolic models, the expression of each gene first needs to be binned into
#high/medium/low/nonexpressed categories (there is also a zero category for genes not present in the dataset/model,
#but we'll deal with that when it inevitably arises). The goal here is to import the read count table that Sambhawa
#sent, see if there are natural breaks, and bin each gene appropriately. For now, we can do it for all samples,
#but we will likely have to split it up by sample to make sample-specific models.

#most of this I'm doing in python, but I need R for graphing.
setwd("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_analyzed_data")

expression_data <- read.table("all_samples_protein_coding_subread_counts.txt", header = TRUE, sep="")
names(expression_data)[1] <- "gene_ID"

#subset data into 500-row chunks
expression_subset1 <- sample_n(expression_data, 500)
expression_subset2 <- sample_n(expression_data, 500)
expression_subset3 <- sample_n(expression_data, 500)
expression_subset4 <- sample_n(expression_data, 500)
expression_subset5 <- sample_n(expression_data, 500)
expression_subset6 <- sample_n(expression_data, 500)
expression_subset7 <- sample_n(expression_data, 500)
expression_subset8 <- sample_n(expression_data, 500)
expression_subset9 <- sample_n(expression_data, 500)
expression_subset10 <- sample_n(expression_data, 500)

#make some graphs.
#will need to facet-wrap by sample, probably. Which means we need to convert from wide to long.
expression_subset1_long <- gather(expression_subset1, sample_id, expression_level, B01:s96, factor_key=TRUE)

#all right, now graph
ggplot(expression_subset1_long, aes(x=gene_ID, y=expression_level)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ sample_id)

#too many facets. Pick one patient ID (s28).
expression_subset1_long_s28 <- expression_subset1_long[ which(expression_subset1_long$sample_id=='s28'), ]

ggplot(expression_subset1_long_s28, aes(x= reorder(gene_ID,-expression_level), y=expression_level)) + 
  geom_bar(stat = "identity") +
  coord_flip()
  
expression_subset1_long_B04 <- expression_subset1_long[ which(expression_subset1_long$sample_id=='B04'), ]

ggplot(expression_subset1_long_B04, aes(x= reorder(gene_ID,-expression_level), y=expression_level)) + 
  geom_bar(stat = "identity") +
  coord_flip()

#######################################
####New plan! Need to actually normalize these data. Will work through the tutorial online elsewhere. 

