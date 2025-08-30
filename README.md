# Variant Calling Pipeline for WES in Rare Disease

### Project Overview

This repository contains a **bioinformatics pipeline** for performing **variant calling** on Whole Exome Sequencing (WES) data from a **neonate affected by MOPD-II (a rare genetic disorder)**. The workflow demonstrates how to process raw sequencing reads, perform quality control, align to a reference genome, call variants, annotate and visualise them and generate summary reports.

This project highlights how WES analysis can uncover pathogenic variants in rare diseases.
Such workflows are critical in clinical genomics, where rapid, reproducible pipelines can directly impact diagnosis and treatment.

The pipeline is designed as a **generalizable workflow** for studying **rare disease genomics**, with MOPD-II used as a case study.

---

### Tools & Dependencies

* **SRA Toolkit** (for FASTQ download)
* **FastQC / MultiQC** (for QC)
* **Chrome** `.deb` package (for viewing multiQC HTML reports) 
* **Trimmomatic** (for trimming)
* **BWA** (for alignment)
* **Samtools and Picard** (for BAM processing)
* **GATK** (for variant calling)
* **VEP** (for variant annotation)
* **IGV** (for visualisation)

---

### Dataset

* **Source:** [SRA](https://www.ncbi.nlm.nih.gov/sra)
* **Accession ID:** `SRR18395025` 
* **Type:** Whole Exome Sequencing (WES) FASTQ files
* **Library:** Illumina NovaSeq 6000
* **Organism:** Human (_Homo sapiens_)
* **Reference Genome:** GRCh38/hg38 (latest assembly)

---

### How to Run

- Clone this repo: `git clone https://github.com/Naila-Srivastava/Variant_Calling-WES-RareDisease.git`
                 cd `Variant_Calling-WES-RareDisease`

```bash
# Example commands

fastq-dump --split-files --gzip -X 100000 SRRxxxxxxx     # Adjust the subset read counts

fastqc SRRxxxxxxxx_*.fastq.gz
multiqc .

trimmomatic PE -threads 4 -phred33 \
SRRxxxxxxxx_1.fastq.gz SRRxxxxxxxx_2.fastq.gz \
SRRxxxxxxxx_1_paired.fq.gz SRRxxxxxxxx_1_unpaired.fq.gz \
SRRxxxxxxxx_2_paired.fq.gz SRRxxxxxxxx_2_unpaired.fq.gz \
ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

bwa mem -t 4 hg38.fa SRRxxxxxxxx_1.fastq SRRxxxxxxxx_2.fastq > SRRxxxxxxxx_aligned.sam

samtools view -Sb SRRxxxxxxxx_aligned.sam > SRRxxxxxxxx_aligned.bam
samtools sort SRRxxxxxxxx_aligned.bam -o SRRxxxxxxxx_sorted.bam

picard MarkDuplicates \
I=SRRxxxxxxxx_sorted.bam \
O=SRRxxxxxxxx_dedup.bam \
M=SRRxxxxxxxx_dedup.metrics.txt

gatk CreateSequenceDictionary -R hg38.fa     # Create a seq dictionary (1 per reference)
gatk HaplotypeCaller \                       # Call variants
    -R hg38.fa \
    -I SRRxxxxxxxx_dedup.bam \
    -O SRRxxxxxxxx_gatk_variants.vcf.gz
```
---

### Methodology

<img width="263" height="921" alt="image" src="https://github.com/user-attachments/assets/769ff303-ab47-48e8-b103-236178f2282a" />

---

### Features

1. End-to-end workflow: from raw reads to annotated VCF.
2. Automated QC reporting using FastQC + MultiQC.
3. Rare disease relevance: pipeline tailored for rare phenotype exploration.
4. Reproducibility: modular scripts for easy reuse across datasets.
5. Organized repo: clear separation of data, scripts, results, and reports.

---

### Visualisations

- FastQC reports → per-sample quality check (per-base quality, GC content, adapter presence).
- MultiQC dashboard → aggregate sample QC metrics in one interactive report.
- Variant summary plots → variant counts, Ti/Tv ratios, functional impact categories.
- Annotated variants table → candidate mutations with gene names, predicted effects, and clinical significance.

---

### Results

* **FastQC + MultiQC reports**
  * Sequencing quality, GC content, adapter content
* **Alignment stats**
  * % mapped reads, coverage
* **Variant files**
  * Raw and filtered VCF files
* **Rare variants identified**
  * Highlight potential pathogenic variants linked to MOPD-II

---

### Key Takeaways

- Rare disease-focused WES pipeline from raw FASTQ → annotated variants.
- Mastery of Linux command-line tools for real-world bioinformatics.
- Clinical insights through variant annotation & prioritization.
- Transferable skills for genomics research & clinical diagnostics.

---

### What's Next

* Expand pipeline to **multiple rare diseases** datasets
* Integrate **Nextflow/Snakemake** for workflow automation
* Add **structural variant detection**
* Perform **functional annotation** and pathway analysis

---

### References

* EBI ENA: https://www.ebi.ac.uk/ena/browser/
* FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* MultiQC: https://multiqc.info/
* Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
* BWA-MEM: http://bio-bwa.sourceforge.net/
* Samtools: https://www.htslib.org/
* Picard: https://broadinstitute.github.io/picard/
* GATK HaplotypeCaller: https://gatk.broadinstitute.org/
* VEP: https://useast.ensembl.org/info/docs/tools/vep/index.html

---

### Project Structure

```plaintext
variant_calling-WES-RareDisease/
│── data/               # FASTQ and reference files (not uploaded due to size)  
│── scripts/            # Bash scripts for each step of the pipeline  
│── results/            # Output files (BAM, VCF, reports)  
│── reports/            # MultiQC & FastQC HTML reports  
│── README.md           # Project documentation (this file)  
