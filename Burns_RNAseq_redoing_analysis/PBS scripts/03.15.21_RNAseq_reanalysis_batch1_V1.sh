#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH --array=1-88%8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

module load kallisto
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B01_S1_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B01_S1_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B01_S1_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B02_S2_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B02_S2_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B02_S2_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B03_S3_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B03_S3_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B03_S3_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B04_S4_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B04_S4_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B04_S4_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B05_S5_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B05_S5_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B05_S5_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B06_S6_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B06_S6_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B06_S6_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B07_S7_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B07_S7_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B07_S7_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B08_S8_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B08_S8_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B08_S8_R2_001.fastq


