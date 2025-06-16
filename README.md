# OxStressLAB-RNAseq: Automated RNA-seq Pipeline for Oxidative Stress Analysis in Lactic Acid Bacteria

**Description:**  
This project implements an automated RNA-seq analysis pipeline using Snakemake to assess the impact of oxidative stress on gene expression in lactic acid bacteria.

---

## Системные требования
| Компонент | Минимальные требования |
|-----------|------------------------|
| ОС | Windows 10+/macOS 10.15+/Linux Ubuntu 20.04+ |
| R | версия 4.2.0 или новее |
| RStudio | 2023.03+ (рекомендуется) |
| Память | 8 ГБ ОЗУ (16+ ГБ для больших датасетов) |
| Место на диске | 5 ГБ свободного пространства |

---

## Project Structure

- **data/**: Raw FASTQ files (paired-end reads)
- **references/**: Reference genome and annotation
- **scripts/**: Analysis scripts (R, Python)
- **Snakefile**: Main pipeline definition for Snakemake
- **results/**: Outputs (tables, plots, reports)
- **logs/**: Log files and reports (MultiQC, Snakemake logs)

---

## Installation and Setup

### 1. Install Miniconda (if not already installed)

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

text

### 2. Clone the Repository

git clone https://github.com/olga-kuta/OxStressLAB-RNAseq-
cd OxStressLAB-RNAseq-

text

### 3. Create and Set Up the Conda Environment

conda create -n RNA-seq
conda activate RNA-seq
conda install -n RNA-seq bioconda::hisat2
conda install -n RNA-seq bioconda::fastqc
conda install -n RNA-seq bioconda::multiqc
conda install -n RNA-seq gffread
conda install -n RNA-seq bioconda::subread
conda install -n RNA-seq bioconda::samtools

text

### 4. Install R Packages

conda activate RNA-seq
R

text

Inside R, run:

install.packages(c("ggplot2", "pheatmap", "ggrepel", "tidyr", "dplyr", "tibble"))
if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "rtracklayer", "clusterProfiler"))
q()

text

---

## Downloading Data

### 1. Download Transcriptomic Data

cd data
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/008/SRR5578738/SRR5578738_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/008/SRR5578738/SRR5578738_2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/009/SRR5578739/SRR5578739_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/009/SRR5578739/SRR5578739_2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/000/SRR5578740/SRR5578740_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/000/SRR5578740/SRR5578740_2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/005/SRR5578735/SRR5578735_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/005/SRR5578735/SRR5578735_2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/006/SRR5578736/SRR5578736_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/006/SRR5578736/SRR5578736_2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/007/SRR5578737/SRR5578737_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR557/007/SRR5578737/SRR5578737_2.fastq.gz
gunzip *.fastq.gz
cd ..

text

### 2. Download Reference Files

cd references
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/617/525/GCF_001617525.2_ASM161752v2/GCF_001617525.2_ASM161752v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/617/525/GCF_001617525.2_ASM161752v2/GCF_001617525.2_ASM161752v2_genomic.gff.gz
gunzip *.gz
cd ..

text

---

## Running the Pipeline

snakemake -s snakemake.py -c <number_of_cores>

text
Replace `<number_of_cores>` with the number of CPU cores to use (e.g., `4`).

---

## Results

- **Quality Control:** MultiQC reports in `results/multiqc/`
- **Differential Expression:** Tables in `results/diffexp/`
- **Visualizations:** Plots (PCA, heatmaps, volcano plots) in `results/plots/`
- **Functional Enrichment:** GO/KEGG results in `results/enrichment/`

---

## Reproducibility

- **Conda environment** ensures all dependencies are tracked.
- **Snakemake workflow** automates the entire analysis.
- **Version control** via Git.

---

## Additional Notes

- **For best results, run all commands from the root of the repository.**
- **If you encounter any issues, please check the logs in `logs/` and consult the Snakemake documentation.**
- **For more information about the data and analysis, see the project repository and linked publications.**

---
