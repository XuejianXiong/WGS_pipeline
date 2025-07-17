#!/bin/bash
# trim_fastp_custom.sh - Trim specified paired-end samples with fastp

set -euo pipefail

SAMPLES=("SRR062634" "SRR062635" "SRR062637")
INPUT_DIR="Data"
OUTPUT_DIR="Result"
THREADS=4

mkdir -p "$OUTPUT_DIR"

for sample in "${SAMPLES[@]}"; do
  fq1="${INPUT_DIR}/${sample}_1.fastq"
  fq2="${INPUT_DIR}/${sample}_2.fastq"

  if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
    echo "ERROR: Missing files for sample $sample. Skipping."
    continue
  fi

  echo "Trimming sample $sample ..."

  fastp \
    -i "$fq1" \
    -I "$fq2" \
    -o "${OUTPUT_DIR}/${sample}_1.trimmed.fastq" \
    -O "${OUTPUT_DIR}/${sample}_2.trimmed.fastq" \
    --detect_adapter_for_pe \
    --thread "$THREADS" \
    --html "${OUTPUT_DIR}/${sample}_fastp_report.html" \
    --json "${OUTPUT_DIR}/${sample}_fastp_report.json" \
    --report_title "fastp report for $sample"

  echo "Finished trimming $sample."
done

echo "All samples processed."
