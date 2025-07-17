# WGS_pipeline

## Download Data:

✅ 1. Download FASTQ files (full WGS) from 1000 genome using SRA Toolkit


✅ 2. Align to the full genome and extract chr22 only

# Index your hg38 or hg19 reference genome
bwa index reference.fasta

# Align each sample to reference
bwa mem -t 8 reference.fasta SRR622461_1.fastq SRR622461_2.fastq | \
    samtools view -b -o SRR622461.bam -

bwa mem -t 8 reference.fasta SRR622462_1.fastq SRR622462_2.fastq | \
    samtools view -b -o SRR622462.bam -

bwa mem -t 8 reference.fasta SRR622463_1.fastq SRR622463_2.fastq | \
    samtools view -b -o SRR622463.bam -

Sort and index each BAM

samtools sort -o SRR622461.sorted.bam SRR622461.bam
samtools index SRR622461.sorted.bam

samtools view -b SRR622461.sorted.bam chr22 > SRR622461.chr22.bam
samtools index SRR622461.chr22.bam

Repeat for 622462 and 622463


##  Quality Contro:

✅ 3. Quality Control (Raw FASTQ)

Tool: FastQC, MultiQC

Output: per-sample QC reports

✅ 4. Trimming (optional)

Tool: Trimmomatic or fastp

Output: cleaned FASTQ files

## Alignment:

✅ 5. Alignment

Tool: BWA-MEM

Output: sorted BAM files

✅ 6. Post-processing

Tools: SAMtools, Picard

Mark duplicates, sort, index, etc.

✅ 7. Base Quality Score Recalibration (BQSR) (optional)

Tool: GATK BaseRecalibrator

Requires known sites: dbSNP, Mills, etc.

## Variant Calling:

✅ 8. Variant Calling

Tool: GATK HaplotypeCaller (SNVs/Indels)

Tool: Manta (SVs)

✅ 9. Variant Filtering

Tool: GATK VariantFiltration, bcftools

## Variant Annotation and Visualization:

✅ 10. Annotation

Tool: Ensembl VEP or SnpEff

✅ 11. Visualization

Tools: IGV, MultiQC, custom Python/R plots


## 💡 Tips for Scaling to AWS

Store raw data & references in S3

Use modular Docker images (as we discussed)

Use WDL/Nextflow for orchestration

Leverage AWS Batch for compute scaling

## Downstream Analysis:

Use GSEA/DESeq2 for downstream if comparing tumor/normal

