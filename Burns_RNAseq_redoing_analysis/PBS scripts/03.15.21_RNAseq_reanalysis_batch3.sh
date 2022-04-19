#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

module load kallisto
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B17_S17_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B17_S17_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B17_S17_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B18_S18_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B18_S18_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B18_S18_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B19_S19_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B19_S19_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B19_S19_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B20_S20_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B20_S20_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B20_S20_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B21_S21_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B21_S21_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B21_S21_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B22_S22_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B22_S22_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B22_S22_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B23_S23_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B23_S23_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B23_S23_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B24_S24_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B24_S24_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B24_S24_R2_001.fastq
