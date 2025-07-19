# 🧬 WGS_pipeline: A Whole Genome Sequencing Analysis Pipeline

This repository provides a modular and reproducible pipeline for analyzing short-read Whole Genome Sequencing (WGS) data. The pipeline processes raw FASTQ files to high-confidence, annotated variants, following widely accepted best practices. It integrates quality control, adapter trimming, alignment, variant calling, and annotation steps using well-established open-source tools.

---

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
- **VEP** – Variant annotation (Ensembl Variant Effect Predictor)  
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
```

---

## 📂 Folder Structure

```
WGS_pipeline/
├── main.wdl                # calls sub-wdls
├── WDL/
│   ├── main_inputs.json      # input json
│   ├── qc_reads.wdl          # fastqc, multiqc, fastp
│   ├── alignment.wdl         # bwa mem + samtools
├── Docker/                   # Dockfiles for based image and other modular images
│   ├── Dockerfile.base         
│   ├── Dockerfile.qc_reads     
├── Data/                     # Raw FASTQ files and reference genome
│   ├── chr22.fa*
│   ├── sample_1.fastq
│   ├── sample_2.fastq
├── Result/                   # Output: trimmed files, BAMs, VCFs ...
├── Report/                   # FastQC and MultiQC reports
├── Scripts/                  # Wrapper scripts for each step
├── env/                      # (Optional) virtual environment
├── README.md                 # Project documentation
```

---

## 🧪 Key Results

After successful execution, the pipeline will generate:

✔️ Trimmed FASTQ files in Result/

✔️ High-quality aligned BAM files

✔️ Raw and filtered VCF files

✔️ Annotated variant reports using VEP

✔️ Quality reports in HTML and JSON via FastQC and MultiQC

We can visually inspect BAMs and VCFs using IGV.

---

## 📘 License

MIT License – feel free to use, adapt, and share.
