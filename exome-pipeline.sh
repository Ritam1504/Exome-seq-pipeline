#!/bin/bash
# Exome Sequencing Pipeline - Ritam Halder
# Usage: bash exome_pipeline.sh

# ========= STEP 0: Setup =========
mkdir -p data results qc trimmed reference

# Download reference genome (GRCh38 small subset for demo)
echo ">>> Downloading reference genome..."
wget -O reference/GRCh38.fa.gz ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip reference/GRCh38.fa.gz

# Index reference for BWA, SAMtools, and GATK
echo ">>> Indexing reference genome..."
bwa index reference/GRCh38.fa
samtools faidx reference/GRCh38.fa
gatk CreateSequenceDictionary -R reference/GRCh38.fa -O reference/GRCh38.dict

# Download demo FASTQ files (small subset from SRA)
echo ">>> Downloading demo FASTQ data..."
wget -O data/sample_R1.fastq.gz https://github.com/biobenkj/sample-data/raw/main/reads/sample_R1.fastq.gz
wget -O data/sample_R2.fastq.gz https://github.com/biobenkj/sample-data/raw/main/reads/sample_R2.fastq.gz

# ========= STEP 1: Quality Control =========
echo ">>> Running FastQC..."
fastqc data/sample_R1.fastq.gz data/sample_R2.fastq.gz -o qc/

echo ">>> Running MultiQC..."
multiqc qc/ -o qc/

# ========= STEP 2: Trimming =========
echo ">>> Trimming adapters with Trimmomatic..."
trimmomatic PE -threads 4 \
  data/sample_R1.fastq.gz data/sample_R2.fastq.gz \
  trimmed/sample_R1.paired.fq.gz trimmed/sample_R1.unpaired.fq.gz \
  trimmed/sample_R2.paired.fq.gz trimmed/sample_R2.unpaired.fq.gz \
  ILLUMINACLIP:/usr/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# ========= STEP 3: Alignment =========
echo ">>> Aligning reads with BWA MEM..."
bwa mem -t 4 reference/GRCh38.fa \
  trimmed/sample_R1.paired.fq.gz trimmed/sample_R2.paired.fq.gz \
  > results/aligned.sam

# ========= STEP 4: Convert, Sort, Index =========
echo ">>> Converting SAM to BAM and sorting..."
samtools view -Sb results/aligned.sam > results/aligned.bam
samtools sort results/aligned.bam -o results/aligned_sorted.bam
samtools index results/aligned_sorted.bam

# ========= STEP 5: Mark Duplicates =========
echo ">>> Marking duplicates with Picard..."
gatk MarkDuplicates \
  -I results/aligned_sorted.bam \
  -O results/aligned_dedup.bam \
  -M results/marked_dup_metrics.txt
samtools index results/aligned_dedup.bam

# ========= STEP 6: Variant Calling =========
echo ">>> Calling variants with GATK HaplotypeCaller..."
gatk HaplotypeCaller \
  -R reference/GRCh38.fa \
  -I results/aligned_dedup.bam \
  -O results/raw_variants.vcf.gz

# ========= STEP 7: Variant Filtering =========
echo ">>> Filtering variants..."
gatk VariantFiltration \
  -R reference/GRCh38.fa \
  -V results/raw_variants.vcf.gz \
  -O results/filtered_variants.vcf.gz \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "BasicFilters"

echo ">>> Pipeline Finished! Outputs are in results/"
