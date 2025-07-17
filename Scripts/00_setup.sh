#!/bin/bash

# ============================================================================
# WGS Pipeline Setup Script
# Author: [Your Name]
# Description: Downloads test FASTQ data, sets up reference genome, and 
#              prepares index files for BWA, Samtools, and GATK.
# Requirements: prefetch, fasterq-dump, bwa, samtools, gatk, curl
# ============================================================================

set -euo pipefail  # Safe Bash settings: exit on error, unset vars are errors, pipe failures are caught

# ------------------------------
# Configurations
# ------------------------------

# Example datasets
# 1000 Genomes Yoruba Trio (only SRR622461 is public)
# SAMPLES=("SRR622461" "SRR622462" "SRR622463")

# 1000 Genomes CEU Trio - too large for testing
# SAMPLES=("SRR622457" "SRR622458" "SRR622459")

# Small test samples from SRA (E. coli)
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

# ------------------------------
# Step 1: Prepare Working Directory
# ------------------------------
mkdir -p "$WORKDIR"
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
  echo "✔ BWA index already exists."
else
  echo "🔧 Creating BWA index..."
  $BWA index "$REF_NAME_FA"
fi

# Samtools FASTA index
if [ -f "${REF_NAME_FA}.fai" ]; then
  echo "✔ Samtools FASTA index already exists."
else
  echo "🔧 Creating FASTA index with Samtools..."
  $SAMTOOLS faidx "$REF_NAME_FA"
fi

# GATK sequence dictionary
DICT_FILE="${REF_NAME_FA%.fa}.dict"
if [ -f "$DICT_FILE" ]; then
  echo "✔ GATK sequence dictionary already exists."
else
  echo "🔧 Creating sequence dictionary with GATK..."
  $GATK CreateSequenceDictionary -R "$REF_NAME_FA"
fi

echo "✅ Reference genome indexing complete."

# ------------------------------
# Done
# ------------------------------
echo "🎉 Setup complete. You can now proceed to alignment!"
