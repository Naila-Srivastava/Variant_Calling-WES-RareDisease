#!/bin/bash

# Step 3: Aligning the reads against the reference genome (pre-download the reference genome)
bwa mem -t 4 hg38.fa SRR18395025_1_paired.fq.gz SRR18395025_2_paired.fq.gz > SRR18395025_aligned.sam
