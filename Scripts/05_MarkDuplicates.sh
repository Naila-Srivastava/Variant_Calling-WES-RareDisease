#!/bin/bash

# Step 5: Mark the duplicates- Picard
picard MarkDuplicates \			     # If Picard is not found, run it as: java -jar picard.jar MarkDuplicates \...
I=SRR18395025_sorted.bam \
O=SRR18395025_dedup.bam \
M=SRR18395025_dedup.metrics.txt

samtools index SRR18395025_dedup.bam
