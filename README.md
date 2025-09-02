# 🧬 WGS_pipeline: A Whole Genome Sequencing Analysis Pipeline

[![Version](https://img.shields.io/badge/version-v1.0.0-blue.svg)](https://github.com/XuejianXiong/WGS_pipeline/releases/tag/v1.0.0)  
[![Python](https://img.shields.io/badge/python-3.13+-brightgreen.svg)](https://www.python.org/)  
[![Docker](https://img.shields.io/badge/docker-latest-blue.svg)](https://www.docker.com/)  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)  

## 📌 Version History

- **[v1.0.0](https://github.com/XuejianXiong/WGS_pipeline/releases/tag/v1.0.0)** – Initial WGS pipeline release using WDL + Docker (**current main branch**)  

---

This repository provides a modular and reproducible pipeline for analyzing short-read Whole Genome Sequencing (WGS) data. The pipeline processes raw FASTQ files to high-confidence variants, following widely accepted best practices. It integrates reads quality control, adapter trimming, alignment, variant calling, and variant quality control steps using well-established open-source tools.


## 📁 Dataset

**Study:**  
The 1000 Genomes Project. *A global reference for human genetic variation*. **Nature** (2015)

**SRA Accessions:**  
- [SRR062634](https://www.ncbi.nlm.nih.gov/sra/SRR062634)
- [SRR062635](https://www.ncbi.nlm.nih.gov/sra/SRR062635)
- [SRR062637](https://www.ncbi.nlm.nih.gov/sra/SRR062637)

**Note:**
- Sample population: Yoruba in Ibadan, Nigeria (YRI)
- Platform: Illumina Genome Analyzer II
- Technology: Paired-end short-read whole-genome sequencing (WGS)
- Objective: Benchmark small variant calling pipeline using high-quality public data
- Reference Genome: A partial reference genome is used to minimize computational load during development and testing. Specifically, chromosome 22 from the hg38 assembly is downloaded from the UCSC Genome Browser.

---

## 🧰 Tech Stack

The pipeline uses a combination of command-line tools and Python-based utilities within a virtual environment:

- **Python**: 3.13.3 (via `venv`)
- **fastqc**: 0.12.1 — for read quality control
- **multiqc**: 1.30 — to aggregate and summarize QC reports
- **fastp**: 1.0.1 — for read trimming and filtering
- **BWA-MEM** – Alignment to the reference genome  
- **SAMtools** – File conversion and sorting  
- **GATK** – Duplicate marking, BQSR, and variant calling  
- **bcftools** – Variant filtering and statistics  
- **IGV (optional)** – Manual visualization of alignments and variants  

---

## 🚀 How to Run the Pipelines

1. **Clone the repository**
```bash
git clone https://github.com/XuejianXiong/WGS_pipeline.git
cd WGS_pipeline
```

2. **Install dependencies**   

- Install bioinformatic tools:
```bash
brew install fastqc fastp bwa samtools bcftools
```

- Manually install GATK and VEP

- Install python packages:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

- Run each step:

```bash
./Scripts/00_setup.sh                             # Download and extract read files (.fastq)
./Scripts/01_qc_reads.sh                          # Run Fastqc and Multiqc
./Scripts/02_trim_fastp.sh                        # Trim read files using fastp
miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_1.json
miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_2.json
miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_3.json
miniwdl run WDL/main_filter_variants.wdl --input WDL/main_inputs_filter_variants.json 
```

---

## 📂 Folder Structure

```
WGS_pipeline/
├── requirements.txt                              # required packages
├── WDL/                                          # WDL files 
│   ├── main_inputs_1.json                        # input json of SRR062634
│   ├── main_inputs_2.json                        # input json of SRR062635
│   ├── main_inputs_3.json                        # input json of SRR062637
│   ├── main_inputs_filter_variants.json          # input json of all trio samples for filtering variants
│   ├── main_variant_calling.wdl                  # the main WDL file from qc_reads to variant_calling
│   ├── main_filter_variants.wdl                  # the main WDL file of joint_genotyping, bcftools_merge, and qc_variant steps
│   ├── qc_reads.wdl                    # fastqc + multiqc
│   ├── trim_fastq.wdl                  # fastp
│   ├── alignment.wdl                   # bwa mem + samtools
│   ├── dedup.wdl                       # gatk SortSam + gatk MarkDuplicates
│   ├── qc_gatk.wdl                     # gatk CollectAlignmentSummaryMetrics + gatk CollectInsertSizeMetrics
│   ├── bqsr.wdl                        # tabix + samtools + gatk BaseRecalibrator + gatk ApplyBQSR
│   ├── variant_calling.wdl             # samtools + gatk HaplotypeCaller
│   ├── joint_genotyping.wdl            # gatk CombineGVCFs + gatk GenotypeGVCFs
│   ├── select_snps.wdl                 # gatk SelectVariants + gatk VariantFiltration
│   ├── select_indels.wdl               # gatk SelectVariants + gatk VariantFiltration
│   ├── bcftools_merge.wdl              # bcftools concat
│   ├── qc_variant.wdl                  # bcftools + plot-vcfstats
├── Docker/                   # Dockfiles for based image and other modular images                  
│   ├── Dockerfile.base       # base image         
│   ├── Dockerfile.qc_reads   # QC sub-image with FastQC and MultiQC
│   ├── Dockerfile.fastp      # sub-image with fastp
│   ├── Dockerfile.alignment  # alignment sub-image with BWA
│   ├── Dockerfile.gatk       # sub-image with GATK and R
│   ├── Dockerfile.bcftools   # sub-image with bcftools
├── Data/                     # Raw FASTQ files, reference genome, and known variants
├── Result/                   # Output: trimmed files, BAMs, VCFs ...
├── Report/                   # FastQC and MultiQC reports
├── Scripts/                  # Wrapper scripts for each step
│   ├── 00_setup.sh           # download sample files, reference genome, and known variants
│   ├── 01_qc_reads.sh        # QC analysis with FastQC and MultiQC
│   ├── 02_trim_fastp.sh      # trim reads with fastp
├── README.md                 # Project documentation
```

---

## 🧪 Key Results

After successful execution, the pipeline will generate:

✔️ Trimmed FASTQ files in Result/

✔️ High-quality aligned BAM files

✔️ Raw and filtered VCF files

✔️ Reads quality reports in HTML and JSON via FastQC and MultiQC

✔️ Variants quality reports in TXT and PDF via bcftools

We can visually inspect BAMs and VCFs using IGV.

---

## 📘 License

MIT License – feel free to use, adapt, and share.
