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

#the goal here is to convert raw RNAseq counts to TPM (transcripts per million).
setwd("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_analyzed_data")


#step 1: read in metadata
# Read the sample information into R
sampleinfo <- read.delim("CRC_metadata_rnaseq.txt")
View(sampleinfo)
sampleinfo

