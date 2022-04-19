#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=16gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adamo010@umn.edu
#SBATCH -p blekhman

module load kallisto
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s58_S77_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s58_S77_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s58_S77_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s59_S78_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s59_S78_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s59_S78_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s60_S79_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s60_S79_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s60_S79_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s61_S80_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s61_S80_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s61_S80_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s63_S81_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s63_S81_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s63_S81_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s73_S82_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s73_S82_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s73_S82_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s74_S83_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s74_S83_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s74_S83_R2_001.fastq
kallisto quant -i ~/Burns_RNAseq_reanalysis/human_transcriptome.idx -b 100 -o ~/Burns_RNAseq_reanalysis/s85_S84_output /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s85_S84_R1_001.fastq /home/blekhman/data_release/umgc/hiseq/160614_D00635_0135_AC95MVANXX/Blekhman_Project_018/s85_S84_R2_001.fastq
