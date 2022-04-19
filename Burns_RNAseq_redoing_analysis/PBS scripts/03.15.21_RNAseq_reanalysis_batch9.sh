#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

module load kallisto
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s43_S65_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s43_S65_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s43_S65_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s44_S66_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s44_S66_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s44_S66_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s47_S69_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s47_S69_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s47_S69_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s48_S70_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s48_S70_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s48_S70_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s51_S73_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s51_S73_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s51_S73_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s52_S74_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s52_S74_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s52_S74_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s56_S75_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s56_S75_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s56_S75_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s57_S76_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s57_S76_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s57_S76_R2_001.fastq
