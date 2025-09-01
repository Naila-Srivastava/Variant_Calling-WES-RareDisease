#!/bin/bash

# Step 6: Variant Calling- GATK HaplotypeCaller (more comprehensive and for clinical/rare diseases)
gatk CreateSequenceDictionary -R hg38.fa     # Create a seq dictionary (1 per reference)
gatk HaplotypeCaller \                       # Call variants
    -R hg38.fa \
    -I SRR18395025_dedup.bam \
    -O SRR18395025_gatk_variants.vcf.gz
