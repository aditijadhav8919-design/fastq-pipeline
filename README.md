# NGS Data Quality Control and Sequence Refinement Pipeline

## 1. Project Introduction

This repository provides a complete automated workflow for the processing of Next-Generation Sequencing (NGS) data. It enables pipelines that can be executed by sequencing platforms often contains technical artifacts, such as adapter sequences, low-quality bases, and sequencing errors that can compromise downstream analysis. This project implements a standardized protocol to evaluate, clean, and validate sequencing data, ensuring it meets the rigorous quality standards required for genomic research.

## 2. Pipeline Workflow
```
Input: FASTQ Files (Paired-end)
         ↓
    FastQC (Raw)
         ↓
Quality Trimming (Cutadapt)
         ↓
    FastQC (Trimmed)
         ↓
  Trimmed FASTQ Files
         ↓
    BWA Alignment
         ↓
    SAM to BAM
         ↓
    Sort & Index
         ↓
  Variant Calling
         ↓
Final Quality Reports
```

## 3. Initial Quality Assessment (FastQC)

The first phase of the pipeline involves a comprehensive analysis of the raw sequencing files using FastQC. This step is crucial for establishing a baseline understanding of the data's health before any modifications are made.

**Key Metrics Evaluated:**
- **Per-base Quality Scores:** Identifies positions where base-quality scores drop, indicating if the sequencing chemistry degraded during the run.
- **Sequence Content Analysis:** Monitoring GC content and base distribution helps detect potential sample contamination or sequencing bias.
- **Adapter Identification:** The tool specifically identifies the presence of synthetic DNA adapters used to attach DNA fragments to the flow cell, which must be removed to prevent alignment errors.

## 4. Sequence Trimming and Adapter Removal (Cutadapt)

Following the initial assessment, the pipeline utilizes Cutadapt to perform precise sequence refinement. This stage is the core processing step of the workflow.

**Operations Performed:**
- **Adapter Trimming:** Removes adapter sequences identified in the previous step, ensuring only the biological DNA remains for analysis.
- **Quality Filtering:** Any reads that fall below a minimum length threshold after trimming are discarded to prevent "noise" in mapping, where short sequences match incorrectly across multiple locations on the genome.

## 5. Post-Trimming Quality Validation (FastQC)

After trimming, a second FastQC analysis validates that the cleaning process was successful and the data now meets the quality control standards for downstream analysis. This round of quality control is essential for the scientific verification of the results.

**Comparative Evaluation:** By comparing the second report to the initial one can report, the effectiveness of the trimming process can be directly observed.

## 6. Read Alignment (BWA)

Cleaned reads are aligned to a reference genome using BWA-MEM algorithm, producing SAM files for further processing.

## 7. SAM to BAM Conversion and Sorting

SAM files are converted to compressed BAM format and sorted by coordinate for efficient downstream analysis.

## 8. Variant Calling (BCFtools)

High-quality alignments are used to identify genetic variants (SNPs and indels) relative to the reference genome.

## 9. Variant Filtering

Variants are filtered based on quality metrics to retain only high-confidence calls for biological interpretation.

## 10. Workflow Organization

The project is structured with a modular design to ensure clarity and reproducibility. By separating the main workflow from individual process definitions, the pipeline remains maintainable and extensible.

## Requirements

**Software Dependencies:**
- Nextflow (≥21.04)
- FastQC (≥0.11.9)
- Cutadapt (≥3.4)
- BWA (≥0.7.17)
- SAMtools (≥1.13)
- BCFtools (≥1.13)

**System Requirements:**
- Linux/Unix operating system
- Minimum 8GB RAM
- Sufficient storage for reference genome and results

## Installation

### Step 1: Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### Step 2: Install dependencies using Conda
```bash
conda create -n ngs-pipeline -c bioconda fastqc cutadapt bwa samtools bcftools
conda activate ngs-pipeline
```

### Step 3: Clone this repository
```bash
git clone https://github.com/aditijadhav8919-design/fastq-pipeline.git
cd fastq-pipeline
```

## Usage

### Prepare Input Data

Place your paired-end FASTQ files in the `data/` directory:
- `sample_R1.fastq.gz` (forward reads)
- `sample_R2.fastq.gz` (reverse reads)

Place reference genome in `reference_genome/` directory:
- `reference.fa`

### Run the Pipeline
```bash
nextflow run workflow.nf
```

### Resume Failed Run
```bash
nextflow run workflow.nf -resume
```

## Output Structure

Results are organized in the `results/` directory:
```
results/
├── 01_fastqc_raw/        # Quality reports before trimming
├── 02_trimmed/           # Trimmed FASTQ files
├── 03_fastqc_trimmed/    # Quality reports after trimming
├── 04_aligned/           # Aligned reads (SAM)
├── 05_bam/               # BAM files
├── 06_sorted/            # Sorted BAM files
├── 07_variants/          # Called variants (VCF)
└── 08_filtered/          # Filtered variants
```

## Project Structure
```
fastq-pipeline/
├── modules/              # Pipeline modules
│   ├── fastqc_raw.nf
│   ├── cutadapt.nf
│   ├── fastqc_trimmed.nf
│   ├── bwa_align.nf
│   ├── sam_to_bam.nf
│   ├── sort_bam.nf
│   ├── index_bam.nf
│   ├── variant_calling.nf
│   └── filtering.nf
├── data/                 # Your input FASTQ files
├── reference_genome/     # Your reference genome
├── results/              # Pipeline outputs
├── workflow.nf           # Main pipeline script
└── nextflow.config       # Configuration
```

## Monitoring Pipeline Progress
```bash
tail -f .nextflow.log
```

## Troubleshooting

**Pipeline fails?**
```bash
nextflow run workflow.nf -resume
```

**Permission error?**
```bash
chmod +x workflow.nf
```

**Out of memory?**
Edit `nextflow.config` and increase memory allocation.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

MIT License

## Contact

GitHub: [@aditijadhav8919-design](https://github.com/aditijadhav8919-design)

## Acknowledgments

- FastQC developers
- Cutadapt developers
- BWA developers
- SAMtools/BCFtools developers
- Nextflow community
