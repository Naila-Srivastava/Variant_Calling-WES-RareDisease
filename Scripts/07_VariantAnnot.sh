#!/bin/bash

# Step 7: Variant annotation- Ensembl VEP
vep -i SRR18395025_gatk_variants.vcf.gz -o SRR18395025_vep_annotated_online.vcf \
  --assembly GRCh38 --everything --database
