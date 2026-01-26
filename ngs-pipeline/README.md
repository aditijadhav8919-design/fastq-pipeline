# NGS Data Quality Control and Sequence Refinement Pipeline

Bioinformatics pipeline for NGS data quality control and adapter trimming using FastQC and Cutadapt

## Overview

This Nextflow pipeline automates the complete workflow for Next-Generation Sequencing (NGS) data analysis, including:
- Raw data quality control (FastQC)
- Adapter trimming (Cutadapt)
- Post-trimming quality assessment
- Read alignment (BWA)
- SAM to BAM conversion
- BAM sorting and indexing
- Variant calling
- Variant filtering

## Directory Structure
```
ngs-pipeline/
├── modules/              # Nextflow process modules
│   ├── fastqc_raw.nf
│   ├── cutadapt.nf
│   ├── fastqc_trimmed.nf
│   ├── bwa_align.nf
│   ├── sam_to_bam.nf
│   ├── sort_bam.nf
│   ├── index_bam.nf
│   ├── variant_calling.nf
│   └── filtering.nf
├── data/                 # Input FASTQ files
├── reference_genome/     # Reference genome files
├── results/              # Pipeline outputs
├── main.nf              # Main pipeline script
├── workflow.nf          # Workflow definition
└── nextflow.config      # Configuration file
```

## Requirements

### Software Dependencies
- Nextflow (≥21.04)
- FastQC (≥0.11.9)
- Cutadapt (≥3.4)
- BWA (≥0.7.17)
- SAMtools (≥1.13)
- BCFtools (≥1.13)

### System Requirements
- Linux/Unix operating system
- Minimum 8GB RAM
- Sufficient storage for reference genome and results

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/ngs-pipeline.git
cd ngs-pipeline
```

### 2. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 3. Install dependencies using Conda (recommended)
```bash
conda create -n ngs-pipeline -c bioconda fastqc cutadapt bwa samtools bcftools
conda activate ngs-pipeline
```

## Usage

### Basic Usage
```bash
nextflow run main.nf
```

### With Custom Parameters
```bash
nextflow run main.nf \
    --reads 'data/*_{R1,R2}.fastq.gz' \
    --reference 'reference_genome/hg38.fa' \
    --outdir 'results'
```

### Resume Failed Run
```bash
nextflow run main.nf -resume
```

## Input Data

Place your paired-end FASTQ files in the `data/` directory with naming convention:
- `sample_R1.fastq.gz` (forward reads)
- `sample_R2.fastq.gz` (reverse reads)

Reference genome should be in `reference_genome/` directory.

## Output Structure
```
results/
├── 01_fastqc_raw/        # Raw data QC reports
├── 02_trimmed/           # Adapter-trimmed reads
├── 03_fastqc_trimmed/    # Post-trimming QC reports
├── 04_aligned/           # Aligned SAM files
├── 05_bam/               # BAM format alignments
├── 06_sorted/            # Sorted BAM files
├── 07_variants/          # Called variants (VCF)
└── 08_filtered/          # Filtered variants
```

## Pipeline Steps

1. **Raw FastQC**: Quality assessment of raw sequencing reads
2. **Adapter Trimming**: Remove adapters and low-quality bases using Cutadapt
3. **Trimmed FastQC**: Quality assessment after trimming
4. **Alignment**: Map reads to reference genome using BWA-MEM
5. **SAM to BAM**: Convert SAM to compressed BAM format
6. **Sorting**: Sort BAM files by coordinates
7. **Indexing**: Create BAM index files
8. **Variant Calling**: Call variants using BCFtools
9. **Filtering**: Filter variants based on quality metrics

## Configuration

Edit `nextflow.config` to customize:
- Resource allocation (CPU, memory)
- Output directories
- Tool-specific parameters

## Monitoring

View pipeline progress:
```bash
tail -f .nextflow.log
```

## Troubleshooting

### Common Issues

**Issue**: Out of memory errors
```bash
# Increase memory in nextflow.config
process.memory = '16 GB'
```

**Issue**: Permission denied
```bash
chmod +x main.nf workflow.nf
```

**Issue**: Resume not working
```bash
# Clean work directory
rm -rf work/
nextflow run main.nf
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

This project is licensed under the MIT License - see LICENSE file for details.

## Citation

If you use this pipeline, please cite:
```
[Your Name]. (2026). NGS Data Quality Control and Sequence Refinement Pipeline. 
GitHub: https://github.com/yourusername/ngs-pipeline
```

## Contact

For questions or issues, please open an issue on GitHub or contact:
- Name: [Your Name]
- Email: [Your Email]

## Acknowledgments

- FastQC developers
- Cutadapt developers
- BWA developers
- SAMtools/BCFtools developers
- Nextflow community
