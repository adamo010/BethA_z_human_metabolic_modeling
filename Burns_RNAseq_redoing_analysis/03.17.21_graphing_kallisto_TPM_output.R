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

#the goal here is to graph TPM (output from kallisto)
setwd("/Users/adamo010/Documents/CORDA_CRC_human_models/Burns_RNAseq_redoing_analysis")

#step 1: read in the data
files_to_read <- list.files(path = "TPM_abundances_from_kallisto/",pattern = ".tsv$",full.names = T)
all_files <- lapply(files_to_read,function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = TRUE)
})
#this gives us a list called all_files, which is a list of all the TPM abundance files

#step 2: log transform (after adding 1, which I... did not know about)
for (file in all_files){
  file$tpm = file$tpm + 1
  #file$log_tpm=log(file$tpm)
  #file[file == -Inf] <- 0
}
#this creates a new variable, log_tpm, in the dataframe file, which is equivalent to log(tpm) in the dataframe file.
#the second line gets rid of the -Inf values that happen when you take the log of 0 TPM; replaces them with 0. Don't need to do this after adding one to all files
#note that, for some reason, the imported dataframes did not conserve capitalization, so TPM became tpm.

#step 3: graph
#as a trial, I'm going to import and graph one dataframe. 
B02_S2 <- read.table(file = 'TPM_abundances_from_kallisto/B02_S2_abundance.tsv', sep = '\t', header = TRUE)
B02_S2$tpm = B02_S2$tpm + 1
B02_S2$log_tpm=log(B02_S2$tpm)
#B02_S2[B02_S2 == -Inf] <- 0

#aha, actually want a histogram: "analyze a histogram of log-transformed TPM values in a sample. "
ggplot(B02_S2, aes(x=log_tpm, fill=cut(log_tpm, 500))) + 
  geom_histogram(bins=250, show.legend = FALSE) +
  ylim(c(0,3000)) +
  scale_x_continuous(breaks=seq(0,10,0.5), limits=c(0,10)) +
  #scale_y_continuous(breaks=seq(0,2000,250), name= "gene counts") +
  #xlim(c(0,10)) +
  ggtitle("B02_S2, binwidth=250") +
  theme(plot.title = element_text(hjust=0.5)) 

#ugh, how do I... save all these graphs as the names of their files, and also change the y axis.
plots <- lapply(all_files,function(x) {
  ggplot(x, aes(x=log_tpm, fill=cut(log_tpm, 500))) + 
    geom_histogram(bins=250, show.legend = FALSE) +
    ylim(c(0,3000)) +
    scale_x_continuous(breaks=seq(0,10,0.5), limits=c(0,10)) +
    ggtitle("B02_S2, binwidth=250") +
    theme(plot.title = element_text(hjust=0.5)) 
})

lapply(names(plots), 
       function(x) ggsave(filename=paste(x,".jpeg",sep=""), plot=plots[[x]]))

for (df in all_files) { 
  histogram <- ggplot(df, aes(x=log_tpm, fill=cut(log_tpm, 500))) + 
    geom_histogram(bins=250, show.legend = FALSE) +
    ylim(c(0,3000)) +
    scale_x_continuous(breaks=seq(0,10,0.5), limits=c(0,10))
  #ggsave(df, plot= histogram, device= "pdf")
}
#ggsave(plots,filename=paste("myplot",nm[i],".png",sep=""))

#this is why I fucking hate R. How hard is it to iterate through a list of dataframes, generate the same fucking
#graph for each dataframe, and save the output? Apparently very fucking hard, and if you don't know how to do it
#already, you're a fuckwit, according to the R community.

pdf("rplots.pdf") 
print(plots)
dev.off() 

###############FUUUUUUUUUUUUUUUUUCK
#step 1: read in the data- have done sample manipulation in python because fukc R.
files_to_read2 <- list.files(path = "TPM_abundances_from_kallisto/",pattern = "adjusted",full.names = T)
all_files2 <- lapply(files_to_read2,function(x) {
  read.table(file = x, 
             sep = ',', 
             header = TRUE)
})

#step 2: plot
plots <- lapply(all_files2,function(x) {
  ggplot(x, aes(x=log_tpm, fill=cut(log_tpm, 500))) + 
    geom_histogram(bins=250, show.legend = FALSE) +
    ylim(c(0,3000)) +
    scale_x_continuous(breaks=seq(0,10,0.5), limits=c(0,10)) +
    ggtitle(x, "binwidth=250") +
    theme(plot.title = element_text(hjust=0.5)) 
})

#step 3: save
pdf("rplots2.pdf") 
print(plots)
dev.off() 

#well... they're all named the same thing, which is unfortunate
