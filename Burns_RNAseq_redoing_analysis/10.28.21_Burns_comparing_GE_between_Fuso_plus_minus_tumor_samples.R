library("ggplot2")
library("DESeq2")
library("dplyr")

#the goal is to import the kallisto-run Burns RNA seq data and run it in DEseq2. The goal output is differentially expressed genes SPECIFICALLY
#between Fusobacterium positive and negative tumor samples. More specifically, get log(fold change) in gene expression here.

#starting point is 04.05.21_DEseq2_analysis_Burns_data_clean.R

#Step 0: clear out environment and set wd
rm(list = ls()) #clear out environment as needed. 
setwd("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis/") 

##############step 1: import counts data (from tximport, or from Sambhawa who has given me the counts table) and metadata
countData <- read.csv('/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_analyzed_data/all_samples_protein_coding_subread_counts.txt', 
                      header = TRUE, sep="\t")

metaData1 <- read.csv("/Users/adamo010/Documents/MICOM_CRC_FBA/Individual_specific_microbiomes/Analyzing_MICOM_communities_Burns2015/10.28.21_key_for_identifying_Fuso_positive_tumor_samples_Burns2015.csv")
#metaData1$Patient_Blind_ID <- as.factor(metaData$Patient_Blind_ID) #set patient_blind_id as factor
metaData1 = subset(metaData1, select= -c(X)) #drop index column X
metaData2 <- read.csv("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_analyzed_data/CRC_metadata_rnaseq.txt", header = TRUE, sep = "\t")
metaData2 <- subset(metaData2, select = c("Tissue.RNA.DNA_Tube_ID", "SampleID")) #only keep relevant columns
metaData <- merge(x = metaData1, y = metaData2, by = "SampleID", all = TRUE) #merge dataframes by common column, SampleID
metaData <- metaData %>%  rename(tube_id = Tissue.RNA.DNA_Tube_ID)

#set tube_id as row names: A Process.
#I don't know why we were using sampleid before, since sampleid doesn't correpsond to anything in countData
#metaData_sub <- metaData[, c("tube_id")] #copy row names as vector
#row.names(metaData) <- metaData_sub #set dataframe row names NOTE: for a dataframe, need to use row.names, not rownames. FUCK R.
#metaData = subset(metaData, select= -c(tube_id)) #remove sample ID as a column, since it's a row name now.

#############step 2: subset countData by tumor samples only.
tumor_metadata <- filter(metaData, Description== "tumor") #subset metadata by tumor only
tumor_tubeids <- pull(tumor_metadata, tube_id) #get a list of tumor tube ids
rownames(countData) <- countData$X #keep gene names as row names
tumor_countData <- countData[ ,which((names(countData) %in% tumor_tubeids)==TRUE)] #only include columns whose names match the tubeids in tumor_tubeids

############step 3: construct deseq2 object:
dds <- DESeqDataSetFromMatrix(countData=tumor_countData, 
                              colData=tumor_metadata, 
                              design = ~Fuso_present) #this is the tricky part. No longer using paired data, so maybe less complicated?

#got a warning message that: In DESeqDataSet(se, design = design, ignoreRank)::some variables in design formula are characters, converting to factors"
#which should be fine
#NOTE here, I"m not filtering at all. Shouldn't be a problem (I hope)

#############step 4: run deseq2 proper:
dds_from_deseq <- DESeq(dds)
#got "5 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest"
#I suspect this is because I did not abundance filter (Samantha had the same problem and the lab agreed this was probably the reason)
#check: 
colData(dds)

#############step 5: look at results:
res <- results(dds_from_deseq)
summary(res) #summary of results for differential gene expression
res <- res[order(res$padj),] #re-order by p-value
head(res)

summary(res) #summary of results for differential gene expression
res <- res[order(res$log2FoldChange),] #re-order by logfoldchange
head(res)

#############step 6: graphing our results:
#plot the top most differentially expressed genes (fill in here)
par(mfrow=c(2,3))
plotCounts(dds_from_deseq, gene="HIST2H4A", intgroup="Fuso_present")
plotCounts(dds_from_deseq, gene="PADI3", intgroup="Fuso_present")
plotCounts(dds_from_deseq, gene="HIST2H3A", intgroup="Fuso_present")
plotCounts(dds_from_deseq, gene="CSAG1", intgroup="Fuso_present")
plotCounts(dds_from_deseq, gene="LCN15", intgroup="Fuso_present")
plotCounts(dds_from_deseq, gene="MKRN3", intgroup="Fuso_present")

#make a volcano plot:
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#make a PCA plot:
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
vsdata <- vst(dds_from_deseq, blind=FALSE)
plotPCA(vsdata, intgroup="Fuso_present") #using the DESEQ2 plotPCA fxn 
plotPCA(vsdata, intgroup=c("Fuso_present", "Patient_Blind_ID"))
#this is neat.
head(vsdata)
head(dds_from_deseq)

############step 6: export the data
#first, subset dataframe to only include significant padj values
resOrdered <- res[order(res$padj),] #re-order by p-value
resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(resSig), file="10.28.21_DEseq2_sig_results_Burns_data_Fuso_presence_absense_tumor_only.csv")

######step 7: vst transform

head(vsdata)
# Convert the DESeq transformed object to a data frame
vsdata <- assay(vsdata) #converts to matrix
vsdata_df <- as.data.frame(vsdata) #converts to dataframe
vsdata_df$Gene <- rownames(vsdata_df) #convert gene names to row names
head(vsdata_df) #cool: each column is a sample_ID, each row is a gene, and it's populated by... ??? hopefully fold change in gene expression?

########step 8: pull out top differentially expressed genes
#pull out the differentially expressed genes and sort by the absolute value of fold change; keep the top 250.
summary(res) #summary of results for differential gene expression
res <- res[order(res$padj),] #re-order by p-value
head(res) #take a look
res <- assay(res) #converts to matrix
res_df <- as.data.frame(res) #converts to dataframe
res_df$Gene <- rownames(res_df) #convert gene names to row names
sigGenes <- rownames(res_df[res_df$padj <= .05,]) #pull out a list of genes (which correspond to rownames) where padj is <=0.05
#gives us a list of 1277 genes.
#pull top 250 most significant genes (smallest p-values)
summary(res) #summary of results for differential gene expression
res <- res[order(res$padj),] #re-order by p-value
head(res) #take a look
res <- assay(res) #converts to matrix
res_df <- as.data.frame(res) #converts to dataframe
top_100_res_df <- res_df %>% top_n(-100, wt= padj)

write.csv(as.data.frame(top_100_res_df), file="10.28.21_DEseq2_sig_results_Burns_data_Fuso_presence_absense_tumor_only_top_100_DE_genes.csv")
