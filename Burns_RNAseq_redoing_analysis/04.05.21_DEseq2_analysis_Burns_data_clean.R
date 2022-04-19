library("ggplot2")
library("DESeq2")
library("dplyr")

#the goal is to import the kallisto-run Burns RNA seq data and run it in DEseq2. The goal output is differentially expressed genes and, 
#more specifically, get log(fold change) in gene expression between cancer and normal tissues. This will then be correlated with metabolic model-based
#microbial growth rates to (hopefully) show something interesting. 

setwd("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_analyzed_data") 

##############step 1: import counts data (from tximport, or from Sambhawa who has given me the counts table) and metadata
countData <- read.csv('all_samples_protein_coding_subread_counts.txt', header = TRUE, sep="\t")
#head(countData)
metaData <- read.csv('CRC_metadata_rnaseq.txt', header = TRUE, sep = "\t")
head(metaData)
#may need to set patient_blind_id as factor
metaData$Patient_Blind_ID <- as.factor(metaData$Patient_Blind_ID)

##ADDED 04.13.21L drop Tissue.RNA.DNA_Tube_ID, since it is unhelpful. Also drop Height, Weight, BMI, RNA_RIN, RNA_conc
metaData = subset(metaData, select = -c(Tissue.RNA.DNA_Tube_ID, Height, Weight, BMI, RNA_RIN, RNA_conc.) )
##ADED 04.16.21: dropping everything in metadata that's not specifically necessary
metaData2 = subset(metaData, select= c(Patient_Blind_ID, SampleID, Description))
class(metaData2) #check the type of object this is; it is a data.frame
row.names(metaData2) #prints row names; very unhelpful

#set SampleID as row names
metaData_sub <- metaData[, c("SampleID")] #copy row names as vector
row.names(metaData) <- metaData_sub #set dataframe row names NOTE: for a dataframe, need to use row.names, not rownames. FUCK R.
metaData = subset(metaData, select= -c(SampleID)) #remove sample ID as a column, since it's a row name now.

metaData2_sub <- metaData2[, c("SampleID")] #copy row names as vector
row.names(metaData2) <- metaData2_sub #set dataframe row names NOTE: for a dataframe, need to use row.names, not rownames. FUCK R.
metaData2 = subset(metaData2, select= -c(SampleID)) #remove sample ID as a column, since it's a row name now.
#continue onto step 2

############step 2: construct deseq2 object:
#have to set first column of count_data as index?
countData_sub <- countData[, c("X")] #copy row names as vector
row.names(countData) <- countData_sub #set dataframe row names NOTE: for a dataframe, need to use row.names, not rownames. FUCK R.
countData2 = subset(countData, select= -c(X)) #remove sample ID as a column, since it's a row name now.

#try again:
dds <- DESeqDataSetFromMatrix(countData=countData2, 
                              colData=metaData, 
                              design = ~Patient_Blind_ID + Description) #this part indicates that we have paired data.

#got a warning message that: In DESeqDataSet(se, design = design, ignoreRank)::some variables in design formula are characters, converting to factors"
#which should be fine
#NOTE here, I"m not filtering at all. Shouldn't be a problem (I hope)

#############step 3: run deseq2 proper:
dds_from_deseq <- DESeq(dds)
#got "5 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest"
#I suspect this is because I did not abundance filter (Samantha had the same problem and the lab agreed this was probably the reason)

#check: 
colData(dds)

#############step 4: look at results:
res <- results(dds_from_deseq)
head(results(dds_from_deseq, tidy=TRUE)) #let's look at the results table

summary(res) #summary of results for differential gene expression
res <- res[order(res$padj),] #re-order by p-value
head(res)

summary(res) #summary of results for differential gene expression
res <- res[order(res$log2FoldChange),] #re-order by p-value
head(res)

#try to see by patient_blind_ID
res_metadata <- mcols(res, use.names=TRUE) #nothing notable

#try to see by patient_blind_ID
res_metadata <- mcols(res, use.names=TRUE) #nothing notable
res_metadata

#############step 5: graphing our results:
#plot the top most differentially expressed genes (fill in here)
par(mfrow=c(2,3))
plotCounts(dds_from_deseq, gene="KRT80", intgroup="Description")
plotCounts(dds_from_deseq, gene="OTOP2", intgroup="Description")
plotCounts(dds_from_deseq, gene="CDH3", intgroup="Description")
plotCounts(dds_from_deseq, gene="BMP3", intgroup="Description")
plotCounts(dds_from_deseq, gene="ATP1A2", intgroup="Description")
plotCounts(dds_from_deseq, gene="AQP8", intgroup="Description")

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
plotPCA(vsdata, intgroup="Description") #using the DESEQ2 plotPCA fxn 
plotPCA(vsdata, intgroup=c("Description", "Patient_Blind_ID"))
#this is neat.
head(vsdata)
head(dds_from_deseq)

############step 6: export the data
#first, subset dataframe to only include significant padj values
resOrdered <- res[order(res$padj),] #re-order by p-value
resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(resSig), file="04.16.21_DEseq2_sig_results_Burns_data_fixed_design.csv")

#Save deseq2 object
save(res, file="04.16.21_DE.RData") #not that I think this will actually save anything, but one can hope

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
#gives us a list of 3600 genes.
#pull top 250 most significant genes (smallest p-values)
summary(res) #summary of results for differential gene expression
res <- res[order(res$padj),] #re-order by p-value
head(res) #take a look
res <- assay(res) #converts to matrix
res_df <- as.data.frame(res) #converts to dataframe
top_250_res_df <- res_df %>% top_n(-250, wt= padj)

write.csv(as.data.frame(top_250_res_df), file="04.23.21_DEseq2_sig_results_Burns_data_to_250_DEGs.csv")


