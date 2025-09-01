#!/bin/bash

# Step 1: Fetch the dataset from SRA Database (Subset of the data)
fastq-dump --split-files --gzip -X 100000 SRR18395025


# Step 2: Run QC and MultiQC
fastqc SRR18395025_*.fastq.gz
multiqc .
