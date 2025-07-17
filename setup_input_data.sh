#!/bin/bash

set -euo pipefail

# Sample configuration
## 1000 Genomes Yoruba Trio Samples. Only SRR622461 is public available
#SAMPLES=("SRR622461" "SRR622462" "SRR622463")

## 1000 Genomes CEU Trio (NA12878/NA12891/NA12892). Too big 93GB
#SAMPLES=("SRR622457" "SRR622458" "SRR622459") 

SAMPLES=("SRR062634" "SRR062635" "SRR062637") 


WORKDIR="Data"

# Tools
PREFETCH="prefetch"
FQ_DUMP="fasterq-dump"


echo "📥 Starting WGS download and FASTQ conversion..."
cd "$WORKDIR"

for SAMPLE in "${SAMPLES[@]}"; do
  echo "🔽 Downloading $SAMPLE..."

  if [ ! -f "$SAMPLE/$SAMPLE.sra" ]; then
    $PREFETCH "$SAMPLE"
  else
    echo "✔ $SAMPLE.sra already exists, skipping download."
  fi

  echo "🔄 Converting $SAMPLE.sra to FASTQ..."
  $FQ_DUMP --split-files "$SAMPLE/$SAMPLE.sra"
done

echo "✅ Download and FASTQ conversion complete."
