# FASTQ Processing Pipeline

A comprehensive Nextflow pipeline for processing paired-end FASTQ files through quality control, alignment, and variant calling.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/aditijadhav8919-design/fastq-pipeline.git
cd fastq-pipeline
```

# Overview
This pipeline converts raw sequencing data into high-confidence variant calls through an eight-step workflow:

1. **Initial QC (FastQC)** – Evaluate raw read quality and detect potential issues.
2. **Read Cleaning (Cutadapt)** – Trim adapters and remove low-quality bases.
3. **QC Verification (FastQC)** – Confirm trimming success and improved quality.
4. **Genome Mapping (BWA-MEM)** – Align reads to the reference genome.
5. **File Conversion (Samtools)** – Convert SAM to space-efficient BAM.
6. **Coordinate Sorting (Samtools)** – Sort reads by genomic coordinates.
7. **Variant Discovery (BCFtools)** – Identify SNPs and indels.
8. **Quality Filtering (BCFtools)** – Retain variants that meet quality thresholds.

# Key Features
- Fully automated from raw reads to filtered variants
- Supports multiple samples with parallel execution
- Quality metrics generated at each step
- Standardized outputs (BAM, VCF)
- Cross-platform compatibility
- Environment management via Conda

# Purpose
This pipeline addresses common sequencing issues by:
- Removing low-quality reads
- Accurately mapping reads to the genome
- Distinguishing true variants from errors
- Delivering ready-to-analyze variant files

Produces ready-to-analyze variant files for downstream research applications
### 3. Activate Environment
```bash
conda activate fastq-pipeline
```

## Usage

### Run the Pipeline
```bash
nextflow run main.nf
```

### With Custom Parameters
```bash
nextflow run main.nf \
  --reads "data/*_R{1,2}.fastq.gz" \
  --ref "reference.fa" \
  --outdir "results"
```

## Pipeline Workflow
```
Input: FASTQ Files (Paired-end)
         ↓
    FastQC (Raw)
         ↓
Quality Trimming (Cutadapt)
         ↓
    FastQC (Trimmed)
         ↓
    BWA Alignment
         ↓
    SAM to BAM
         ↓
    Sort & Index
         ↓
  Mark Duplicates
         ↓
  Variant Calling
         ↓
  Filter Variants
         ↓
Final Quality Reports
```

## Requirements

- Conda or Mamba
- At least 8 GB RAM
- Linux or macOS

## Output

Results are saved in the `results/` directory with the following structure:

- `fastqc_raw/` - Raw read QC reports
- `cutadapt/` - Trimmed reads
- `fastqc_trim/` - Trimmed read QC reports
- `alignment/` - BAM files
- `variants/` - VCF files
- `multiqc/` - Aggregated QC report
