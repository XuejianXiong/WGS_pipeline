#!/bin/bash

set -euo pipefail

# Sample IDs
SAMPLES=("SRR062634" "SRR062635" "SRR062637")

# Directory where your FASTQ files are located
FASTQ_DIR="Data"

# Output directory for FastQC and MultiQC results
QC_DIR="Report"

# Create output directory
if [ ! -d "$QC_DIR" ]; then
  echo "📁 Creating QC output directory: $QC_DIR"
  mkdir "$QC_DIR"
else
  echo "ℹ️ QC output directory already exists: $QC_DIR"
fi

echo "🚀 Running FastQC..."

# Run FastQC on all paired FASTQ files
for SAMPLE in "${SAMPLES[@]}"; do
  echo "🔬 Processing $SAMPLE"
  fastqc -o "$QC_DIR" "$FASTQ_DIR/${SAMPLE}_1.fastq" "$FASTQ_DIR/${SAMPLE}_2.fastq"
done

echo "📊 Running MultiQC to summarize FastQC results..."
multiqc "$QC_DIR" -o "$QC_DIR"

echo "✅ QC complete. MultiQC report saved to: $QC_DIR/multiqc_report.html"
