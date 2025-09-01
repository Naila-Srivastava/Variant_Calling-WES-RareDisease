#!/bin/bash

# Step 4: BAM Processing- Samtools
samtools view -Sb SRR18395025_aligned.sam > SRR18395025_aligned.bam	# Convert SAM to BAM
samtools sort SRR18395025_aligned.bam -o SRR18395025_sorted.bam		# Sort the BAM file