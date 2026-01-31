# FASTQ Processing Pipeline

A comprehensive Nextflow pipeline for processing paired-end FASTQ files through quality control, alignment, and variant calling.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/aditijadhav8919-design/fastq-pipeline.git
cd fastq-pipeline
```

### 2. Create Conda Environment
```bash
conda env create -f environment.yml
```

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
