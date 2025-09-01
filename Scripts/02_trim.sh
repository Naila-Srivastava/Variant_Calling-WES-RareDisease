#!/bin/bash

# Step 2: Trimmomatic- To trim adapters & low-base sequences
trimmomatic PE -threads 4 -phred33 \
SRR18395025_1.fastq.gz SRR18395025_2.fastq.gz \
SRR18395025_1_paired.fq.gz SRR18395025_1_unpaired.fq.gz \
SRR18395025_2_paired.fq.gz SRR18395025_2_unpaired.fq.gz \
ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
