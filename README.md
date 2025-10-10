# 🧬 WGS_pipeline: A Whole Genome Sequencing Analysis Pipeline

[![Nextflow](https://img.shields.io/badge/Nextflow-Workflow-orange?logo=nextflow&logoColor=white)](https://www.nextflow.io/)
[![WDL](https://img.shields.io/badge/WDL-Workflow-blue?logo=workflow&logoColor=white)](https://openwdl.org/)  
[![Docker](https://img.shields.io/badge/Docker-Container-blue?logo=docker&logoColor=white)](https://www.docker.com/)  
[![GATK](https://img.shields.io/badge/GATK-Genome%20Analysis-green?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAFElEQVQYV2NkYGD4z0AEMIIBwzAAANcYAR9gCqbRAAAAAElFTkSuQmCC)](https://gatk.broadinstitute.org/)  
[![BWA-MEM](https://img.shields.io/badge/BWA--MEM-Aligner-lightgrey)](http://bio-bwa.sourceforge.net/)  
[![FastQC](https://img.shields.io/badge/FastQC-QC%20Tool-yellow)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)


## 📌 Version History

- **[v2.0.0]** – Refactored WGS pipeline using Nextflow DSL2 + Docker (**current feature/nextflow branch**)  
- **[v1.0.0](https://github.com/XuejianXiong/WGS_pipeline/releases/tag/v1.0.0)** – Initial release: WGS pipeline using WDL + Docker

---

This repository provides a modular and reproducible pipeline for analyzing short-read Whole Genome Sequencing (WGS) data. The pipeline processes raw FASTQ files to high-confidence variants, following widely accepted best practices. It integrates reads quality control, adapter trimming, alignment, variant calling, and variant quality control steps using well-established open-source tools.

The workflow can be executed using either WDL or Nextflow, making it portable across local machines, HPC, or cloud environments.

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

- Manually install GATK and VEP from Broad and Ensembl.

- Install python packages:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

- Run with WDL:

```bash
./Scripts/00_setup.sh                             # Download and extract read files (.fastq)
./Scripts/01_qc_reads.sh                          # Run Fastqc and Multiqc
./Scripts/02_trim_fastp.sh                        # Trim read files using fastp
miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_1.json
miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_2.json
miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_3.json
miniwdl run WDL/main_filter_variants.wdl --input WDL/main_inputs_filter_variants.json 
```

Or run with Nextflow:

```bash
./Scripts/00_setup.sh                             # Download and extract read files (.fastq)
nextflow config                                   # Configure nextflow
nextflow run nextflow/main.nf -params-file nextflow/input.json -profile docker -c nextflow.config
```

---

## 📂 Folder Structure

```
WGS_pipeline/
├── requirements.txt          # Python dependencies
├── .nextflow.config          # Nextflow configuration
├── nextflow/                 # Nextflow DSL2 modules and workflows
├── WDL/                      # WDL workflows and input files 
├── Docker/                   # Dockfiles for based image and other modular images                  
├── Data/                     # Raw FASTQ files, reference genome, and known variants
├── Result/                   # Pipeline outputs (FASTQ, BAM, VCF…)
├── Report/                   # FastQC and MultiQC reports
├── Scripts/                  # Wrapper scripts for each step
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
