#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

module load kallisto
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B09_S9_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B09_S9_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B09_S9_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B10_S10_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B10_S10_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B10_S10_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B11_S11_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B11_S11_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B11_S11_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B12_S12_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B12_S12_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B12_S12_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B13_S13_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B13_S13_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B13_S13_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B14_S14_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B14_S14_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B14_S14_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B15_S15_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B15_S15_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B15_S15_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/B16_S16_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B16_S16_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/B16_S16_R2_001.fastq
