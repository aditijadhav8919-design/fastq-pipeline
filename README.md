# FASTQ Processing Pipeline

A comprehensive Nextflow pipeline for processing paired-end FASTQ files through quality control, alignment, and variant calling.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/aditijadhav8919-design/fastq-pipeline.git
cd fastq-pipeline
```

Overview

This pipeline processes raw sequencing data to generate high-confidence variant calls through an eight-step workflow:

Initial QC (FastQC) – Assess raw read quality and detect potential sequencing issues.

Read Cleaning (Cutadapt) – Trim adapter sequences and remove low-quality bases from read ends.

QC Verification (FastQC) – Confirm trimming success and evaluate improved read quality.

Genome Mapping (BWA-MEM) – Align cleaned reads to the reference genome.

File Conversion (Samtools) – Convert alignment files from SAM to compressed BAM format.

Coordinate Sorting (Samtools) – Sort aligned reads by genomic position.

Variant Discovery (BCFtools) – Identify SNPs and indels in the mapped reads.

Quality Filtering (BCFtools) – Retain variants that meet defined quality thresholds.

Key Features

Fully automated workflow from raw reads to filtered variants

Supports multiple samples with parallel execution

Generates quality metrics at every stage

Produces standardized outputs (BAM and VCF files)

Compatible across platforms

Environment managed via Conda

Purpose:

Sequencing data often contains errors and technical artifacts. This pipeline:

Detects and removes low-quality sequences

Maps reads accurately to their genomic locations

Differentiates true genetic variants from sequencing errors

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
