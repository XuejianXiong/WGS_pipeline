# 🧬 WGS_pipeline: A Whole Genome Sequencing Analysis Pipeline

[![Version](https://img.shields.io/badge/version-2.0-blue.svg)](https://github.com/XuejianXiong/WGS_pipeline/releases)
[![Nextflow](https://img.shields.io/badge/Nextflow-Workflow-orange?logo=nextflow&logoColor=white)](https://www.nextflow.io/)
[![WDL](https://img.shields.io/badge/WDL-Workflow-blue?logo=workflow&logoColor=white)](https://openwdl.org/)  
[![Docker](https://img.shields.io/badge/Docker-Container-blue?logo=docker&logoColor=white)](https://www.docker.com/)  
[![GATK](https://img.shields.io/badge/GATK-Genome%20Analysis-green?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAAFElEQVQYV2NkYGD4z0AEMIIBwzAAANcYAR9gCqbRAAAAAElFTkSuQmCC)](https://gatk.broadinstitute.org/)  
[![BWA-MEM](https://img.shields.io/badge/BWA--MEM-Aligner-lightgrey)](http://bio-bwa.sourceforge.net/)  
[![FastQC](https://img.shields.io/badge/FastQC-QC%20Tool-yellow)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)


## 📌 Version History

- **[v2.0](https://github.com/XuejianXiong/WGS_pipeline)** – Refactored WGS pipeline using **Nextflow DSL2** with full **Docker** support (**current main branch**)
- **[v1.0](https://github.com/XuejianXiong/WGS_pipeline/tree/v1.0)** – Initial release using **WDL + Docker**

---

## 🧠 Overview

**WGS_pipeline** provides a modular, scalable, and fully reproducible workflow for analyzing short-read Whole Genome Sequencing (WGS) data.  
It follows GATK Best Practices to process raw FASTQ reads into high-confidence small variants (SNPs and INDELs).

The pipeline integrates widely used open-source bioinformatics tools and supports execution via either **WDL** or **Nextflow**, enabling flexible deployment on local machines, HPC clusters, or cloud environments.

---

## 📁 Dataset

**Study:**  
The 1000 Genomes Project. *A global reference for human genetic variation*. **Nature** (2015)

**SRA Accessions:**  
- [SRR062634](https://www.ncbi.nlm.nih.gov/sra/SRR062634)
- [SRR062635](https://www.ncbi.nlm.nih.gov/sra/SRR062635)
- [SRR062637](https://www.ncbi.nlm.nih.gov/sra/SRR062637)

**Metadata:**
- **Population:** Yoruba in Ibadan, Nigeria (YRI)  
- **Platform:** Illumina Genome Analyzer II  
- **Technology:** Paired-end short-read WGS  
- **Objective:** Benchmark a small variant calling pipeline using open public data  
- **Reference Genome:** Partial hg38 (chromosome 22) from the UCSC Genome Browser — used to reduce computational load for development and testing

---

## 🧰 Tech Stack

**Core Tools and Languages:**

| Component | Version / Description |
|------------|------------------------|
| **Python** | 3.13.3 (via `venv`) |
| **fastqc** | 0.12.1 — Read quality control |
| **multiqc** | 1.30 — Aggregate QC reports |
| **fastp** | 1.0.1 — Adapter trimming and quality filtering |
| **BWA-MEM** | Short-read alignment |
| **SAMtools** | BAM file processing and sorting |
| **GATK** | Duplicate marking, BQSR, variant calling, and filtering |
| **bcftools** | Variant filtering, merging, and statistics |
| **IGV (optional)** | Visualization of BAMs and VCFs |

---

## 🚀 How to Run the Pipelines

### 1. Clone the Repository

```bash
git clone https://github.com/XuejianXiong/WGS_pipeline.git
cd WGS_pipeline
```

### 2. Build or Pull Docker Images

The pipeline relies on modular Docker containers (e.g., base, qc, bwa, gatk, etc.) for reproducibility and consistent environments.
You can either build them locally or pull pre-built ones if available.


### 2. Install Dependencies

- Bioinformatics tools:

```bash
brew install fastqc fastp bwa samtools bcftools
```

- Install GATK and VEP manually from the Broad Institute and Ensembl, respectively..

- Install python packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### 3. Run with WDL:

```bash
./Scripts/00_setup.sh                             # Download and extract read files (.fastq)
./Scripts/01_qc_reads.sh                          # Run Fastqc and Multiqc
./Scripts/02_trim_fastp.sh                        # Trim read files using fastp

miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_1.json
miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_2.json
miniwdl run WDL/main_variant_calling.wdl --input WDL/main_inputs_3.json
miniwdl run WDL/main_filter_variants.wdl --input WDL/main_inputs_filter_variants.json 
```

### 4. Run with Nextflow:

```bash
./Scripts/00_setup.sh                             # Download and extract read files (.fastq)
nextflow config                                   # Inspect configuration
nextflow run nextflow/main.nf \
          -params-file nextflow/input.json \
          -profile docker \
          -c nextflow.config
```

---

## 📂 Folder Structure

```
WGS_pipeline/
├── requirements.txt          # Python dependencies
├── .nextflow.config          # Nextflow configuration
├── nextflow/                 # Nextflow DSL2 modules and workflows
├── WDL/                      # WDL workflows and input files 
├── Docker/                   # Dockfiles for modular images                 
├── Data/                     # Raw FASTQ files, reference genome, and known variants
├── Result/                   # Pipeline outputs (FASTQ, BAM, VCF…)
├── Report/                   # QC and MultiQC reports
├── Scripts/                  # Wrapper scripts for pipeline steps
├── README.md                 # Project documentation
```

---

## 🧪 Key Results

After successful execution, the pipeline will generate:

✔️ Trimmed FASTQ files

✔️ High-quality aligned BAM files

✔️ Raw and filtered VCFs (SNPs and INDELs)

✔️ Reads QC reports in HTML and JSON from FastQC and MultiQC

✔️ Variants QC statistics and plots from bcftools

We can visually inspect BAMs and VCFs using IGV.

---

## 📘 License

MIT License – feel free to use, adapt, and share.

