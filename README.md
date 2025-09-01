# Variant Calling Pipeline for WES in Rare Disease

_Rare disease genomics with WES: from raw FASTQ to annotated variants_

### Project Overview

An end-to-end **bioinformatics variant calling pipeline** for **Whole Exome Sequencing (WES) data** in **rare disease genomics**, demonstrated using a **MOPD-II neonate** case study. The WES analysis uncovers pathogenic variants in targeted rare diseases.
The workflow demonstrates how to process raw sequencing reads, perform quality control, align to a reference genome, call variants, annotate and visualise them and generate summary reports.
Such workflows are critical in clinical genomics, where rapid, reproducible pipelines can directly impact diagnosis and treatment.

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

This dataset (SRR18395025) from the Sri Lanka Rare Disease Project provides Whole Exome Sequencing of a neonate affected by MOPD-II, a rare genetic disorder with clinical relevance in growth and developmental pathways

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

#### Example commands

```bash
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
gatk HaplotypeCaller \                        # Call variants
    -R hg38.fa \
    -I SRRxxxxxxxx_dedup.bam \
    -O SRRxxxxxxxx_gatk_variants.vcf.gz

vep -i SRRxxxxxxxx_gatk_variants.vcf.gz -o SRRxxxxxxxx_vep_annotated_online.vcf \
  --assembly GRCh38 --everything --database
```
---

### Methodology

<img width="263" height="921" alt="image" src="https://github.com/user-attachments/assets/769ff303-ab47-48e8-b103-236178f2282a" />

---

### Features

1. **End-to-end workflow:** From raw reads to annotated VCF.
2. Automated QC reporting using FastQC + MultiQC.
3. **Rare disease relevance:** Pipeline tailored for rare phenotype exploration.
4. **Reproducibility:** Modular scripts for easy reuse across datasets.
5. **Organised repo:** Clear separation of data, scripts, results, and reports.

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

<img width="1401" height="270" alt="image" src="https://github.com/user-attachments/assets/00c6f110-8f82-4fbb-8f3d-28c97644253d" />

<img width="1398" height="653" alt="image" src="https://github.com/user-attachments/assets/8c476196-3cc9-4efc-987c-a85de949253b" />

<img width="1390" height="572" alt="image" src="https://github.com/user-attachments/assets/b86c1c36-700f-4315-b5b5-432d9683893d" />

<img width="1386" height="607" alt="image" src="https://github.com/user-attachments/assets/8e73f395-6384-4890-95ad-35aca1aed240" />

* **Rare variants identified**: From a total of ~16,250 variants annotated with Ensembl VEP, filtering based on predicted impact (missense, splice-site, nonsense) and computational pathogenicity predictions (SIFT and PolyPhen) resulted in a shortlist of candidate variants of interest.
  
  * Most retained variants are missense mutations with predicted moderate impact.
  * Several variants were consistently annotated as deleterious (SIFT) and probably damaging (PolyPhen), highlighting their potential functional significance.
  * Key affected genes include **MMP23B, CFAP74, PRDM16, CHD5, RNF207, PLEKHG5, DNAJC11, EXOSC10, C1orf167, CLCN6, NPPA, VPS13D**, among others.
  * No strong ClinVar pathogenic annotations were detected in this subset (which might indicate that variants may still be novel).

Summary:
The integration of functional consequence, in silico prediction, and allele effect prioritization identified a set of candidate variants likely to have biological relevance and warrant further investigation.

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

* [EBI ENA](https://www.ebi.ac.uk/ena/browser/)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [BWA-MEM](http://bio-bwa.sourceforge.net/)
* [Samtools](https://www.htslib.org/)
* [Picard](https://broadinstitute.github.io/picard/)
* [GATK HaplotypeCaller](https://gatk.broadinstitute.org/)
* [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)

---

### Project Structure

```plaintext
variant_calling-WES-RareDisease/
│── data/               # FASTQ and reference files (not uploaded due to size)  
│── scripts/            # Bash scripts for each step of the pipeline  
│── results/            # Output files (BAM, VCF, reports)  
│── reports/            # MultiQC & FastQC HTML reports  
│── README.md           # Project documentation (this file)  
