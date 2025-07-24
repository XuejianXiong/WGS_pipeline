#!/bin/bash

# ============================================================================
# WGS Pipeline Setup Script
# Author: Xuejian Xiong
# Description: Downloads test FASTQ data, sets up reference genome, and 
#              prepares index files for BWA, Samtools, and GATK.
# Requirements: prefetch, fasterq-dump, bwa, samtools, gatk, bcftools, 
#               bgzip, tabix, curl, wget
# ============================================================================

# Safe Bash settings: exit on error, unset vars are errors, pipe failures are caught
set -euo pipefail  

# ------------------------------
# Configurations
# ------------------------------

# Example datasets
# 1000 Genomes Yoruba Trio (only SRR622461 is public)
# SAMPLES=("SRR622461" "SRR622462" "SRR622463")

# 1000 Genomes CEU Trio - too large for testing
# SAMPLES=("SRR622457" "SRR622458" "SRR622459")

# Small test samples from SRA
SAMPLES=("SRR062634" "SRR062635" "SRR062637")

# Working directory for input/output
WORKDIR="Data"

# Reference genome (chromosome 22 from hg38) from UCSC
REF_NAME_FA="chr22.fa"
REF_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/${REF_NAME_FA}.gz"

# Tool commands (assumes they're on PATH or configured properly)
PREFETCH="prefetch"
FQ_DUMP="fasterq-dump"
BWA="bwa"
SAMTOOLS="samtools"
GATK="gatk"
BCFTOOLS="bcftools"

# ------------------------------
# Step 1: Prepare Working Directory
# ------------------------------
if [ -d "$WORKDIR" ]; then
  echo "📁 Directory $WORKDIR already exists."
else
  echo "📁 Creating directory $WORKDIR..."
  mkdir -p "$WORKDIR"
fi
cd "$WORKDIR"

# ------------------------------
# Step 2: Download SRA Files & Convert to FASTQ
# ------------------------------
echo "📥 Starting SRA download and FASTQ conversion..."

for SAMPLE in "${SAMPLES[@]}"; do
  echo "🔽 Processing $SAMPLE..."

  # Download .sra file if not already present
  if [ ! -f "$SAMPLE/$SAMPLE.sra" ]; then
    echo "⬇ Downloading $SAMPLE.sra..."
    $PREFETCH "$SAMPLE"
  else
    echo "✔ $SAMPLE.sra already exists."
  fi

  # Convert to paired FASTQ if not already present
  if [[ -f "${SAMPLE}_1.fastq" && -f "${SAMPLE}_2.fastq" ]]; then
    echo "✔ ${SAMPLE}_1.fastq and ${SAMPLE}_2.fastq already exist. Skipping conversion."
  else
    echo "🔄 Converting $SAMPLE.sra to FASTQ..."
    $FQ_DUMP --split-files "$SAMPLE/$SAMPLE.sra"
  fi
done

echo "✅ All samples processed."

# ------------------------------
# Step 3: Download Reference Genome (chr22)
# ------------------------------
echo "📥 Checking reference genome: $REF_NAME_FA..."

if [ ! -f "$REF_NAME_FA" ]; then
  echo "⬇ Downloading $REF_NAME.fa.gz from UCSC..."
  curl -O "$REF_URL"
  gunzip "${REF_NAME_FA}.gz"
else
  echo "✔ $REF_NAME_FA already exists."
fi

# ------------------------------
# Step 4: Index Reference Genome
# ------------------------------

# BWA index
if [ -f "${REF_NAME_FA}.bwt" ]; then
  echo "✔ BWA index of reference already exists."
else
  echo "🔧 Creating reference BWA index..."
  $BWA index "$REF_NAME_FA"
fi

# Samtools FASTA index
if [ -f "${REF_NAME_FA}.fai" ]; then
  echo "✔ Samtools FASTA index of reference already exists."
else
  echo "🔧 Creating reference FASTA index with Samtools..."
  $SAMTOOLS faidx "$REF_NAME_FA"
fi

# GATK sequence dictionary
DICT_FILE="${REF_NAME_FA%.fa}.dict"
if [ -f "$DICT_FILE" ]; then
  echo "✔ GATK sequence dictionary of reference already exists."
else
  echo "🔧 Creating reference sequence dictionary with GATK..."
  $GATK CreateSequenceDictionary -R "$REF_NAME_FA"
fi

echo "✅ Reference genome indexing complete."


# ------------------------------
# Step 5: Download Known Sites for BQSR
# ------------------------------

echo "📥 Preparing known variant sites for BaseRecalibrator..."

# dbSNP for GRCh38 (GATK version)
DBSNP_VCF="dbsnp.vcf"
DBSNP_IDX="dbsnp.vcf.idx"
if [[ -f "$DBSNP_VCF" && -f "$DBSNP_IDX" ]]; then
  echo "✔ dbSNP already downloaded."
else
  echo "⬇ Downloading dbSNP..."
  wget -O dbsnp.vcf https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
  wget -O dbsnp.vcf.idx https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
fi

DBSNP_BGZ="dbsnp.vcf.gz"
DBSNP_TBI="dbsnp.vcf.gz.tbi"
if [[ -f "$DBSNP_BGZ" && -f "$DBSNP_TBI" ]]; then
  echo "✔ dbSNP VCF already compressed and indexed."
else
  echo "🔧 Compressing and indexing dbSNP..."
  bgzip -c "$DBSNP_VCF" > "$DBSNP_BGZ"
  tabix -p vcf "$DBSNP_BGZ"
fi

# Mills and 1000G gold standard indels for GRCh38
MILLS_BGZ="mills.vcf.gz"
MILLS_TBI="mills.vcf.gz.tbi"
if [[ -f "$MILLS_BGZ" && -f "$MILLS_TBI" ]]; then
  echo "✔ Mills and 1000G indels already downloaded."
else
  echo "⬇ Downloading Mills and 1000G indels..."
  wget -O mills.vcf.gz https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
  wget -O mills.vcf.gz.tbi https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
fi

# Subset to chr22
DBSNP_CHR22="dbsnp_chr22.vcf.gz"
MILLS_CHR22="mills_chr22.vcf.gz"

if [[ -f "$DBSNP_CHR22" && -f "$MILLS_CHR22" ]]; then
  echo "✔ Subset VCFs for chr22 already exist."
else
  echo "🧬 Subsetting dbSNP and Mills to chr22..."
  $BCFTOOLS view -r chr22 "$DBSNP_BGZ" -Oz -o "$DBSNP_CHR22"
  tabix -p vcf "$DBSNP_CHR22"

  $BCFTOOLS view -r chr22 "$MILLS_BGZ" -Oz -o "$MILLS_CHR22"
  tabix -p vcf "$MILLS_CHR22"
fi

echo "✅ Known sites for BQSR prepared: $DBSNP_CHR22 and $MILLS_CHR22"


# ------------------------------
# Done
# ------------------------------
echo "🎉 Setup complete. You can now proceed to alignment!"
