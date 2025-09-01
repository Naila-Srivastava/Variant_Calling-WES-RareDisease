#!/bin/bash

# Step 1: Fetch the dataset from SRA Database (Subset of the data)
fastq-dump --split-files --gzip -X 100000 SRR18395025


# Step 2: Run QC and MultiQC
fastqc SRR18395025_*.fastq.gz
multiqc .


# Step 3: Trimmomatic- To trim adapters & low-base sequences
trimmomatic PE -threads 4 -phred33 \
SRR18395025_1.fastq.gz SRR18395025_2.fastq.gz \
SRR18395025_1_paired.fq.gz SRR18395025_1_unpaired.fq.gz \
SRR18395025_2_paired.fq.gz SRR18395025_2_unpaired.fq.gz \
ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# Step 4: Run QC and MultiQC again (optional but recommended)
fastqc SRR18395025_1.fastq.gz SRR18395025_2.fastq.gz \
       SRR18395025_1_paired.fq.gz SRR18395025_2_paired.fq.gz    
multiqc .


# Step 5: Aligning the reads against the reference genome (pre-download the reference genome)
bwa mem -t 4 hg38.fa SRR18395025_1_paired.fq.gz SRR18395025_2_paired.fq.gz > SRR18395025_aligned.sam


# Step 6: BAM Processing- Samtools
samtools view -Sb SRR18395025_aligned.sam > SRR18395025_aligned.bam	# Convert SAM to BAM
samtools sort SRR18395025_aligned.bam -o SRR18395025_sorted.bam		# Sort the BAM file 


# Step 7: Mark the duplicates- Picard
picard MarkDuplicates \			     # If Picard is not found, run it as: java -jar picard.jar MarkDuplicates \...
I=SRR18395025_sorted.bam \
O=SRR18395025_dedup.bam \
M=SRR18395025_dedup.metrics.txt

samtools index SRR18395025_dedup.bam


# Step 8: Variant Calling- GATK HaplotypeCaller (more comprehensive and for clinical/rare diseases)
gatk CreateSequenceDictionary -R hg38.fa     # Create a seq dictionary (1 per reference)
gatk HaplotypeCaller \                       # Call variants
    -R hg38.fa \
    -I SRR18395025_dedup.bam \
    -O SRR18395025_gatk_variants.vcf.gz


# Step 9: Variant annotation- Ensembl VEP
vep -i SRR18395025_gatk_variants.vcf.gz -o SRR18395025_vep_annotated_online.vcf \
  --assembly GRCh38 --everything --database

