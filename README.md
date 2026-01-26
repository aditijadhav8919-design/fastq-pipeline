# NGS Data Quality Control Pipeline

A Nextflow pipeline for NGS data quality control and variant calling.

## What This Pipeline Does

1. **Quality Control** - Checks quality of raw sequencing data (FastQC)
2. **Adapter Trimming** - Removes adapters and low-quality bases (Cutadapt)
3. **Alignment** - Maps reads to reference genome (BWA)
4. **Variant Calling** - Identifies genetic variants (BCFtools)
5. **Filtering** - Filters variants by quality

## Requirements

- Nextflow
- FastQC
- Cutadapt
- BWA
- SAMtools
- BCFtools

## Installation

### Step 1: Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### Step 2: Install dependencies with Conda
```bash
conda create -n ngs-pipeline -c bioconda fastqc cutadapt bwa samtools bcftools
conda activate ngs-pipeline
```

### Step 3: Clone this repository
```bash
git clone https://github.com/aditijadhav8919-design/fastq-pipeline.git
cd fastq-pipeline
```

## How to Run

### Step 1: Prepare your data
Place your FASTQ files in the `data/` folder:
- `sample_R1.fastq.gz` (forward reads)
- `sample_R2.fastq.gz` (reverse reads)

### Step 2: Add reference genome
Place your reference genome in `reference_genome/` folder:
- `reference.fa`

### Step 3: Run the pipeline
```bash
nextflow run workflow.nf
```

## Output

Results will be saved in the `results/` folder:
- `01_fastqc_raw/` - Quality reports before trimming
- `02_trimmed/` - Trimmed FASTQ files
- `03_fastqc_trimmed/` - Quality reports after trimming
- `04_aligned/` - Aligned reads
- `05_bam/` - BAM files
- `06_sorted/` - Sorted BAM files
- `07_variants/` - Called variants
- `08_filtered/` - Filtered variants

## Project Structure
```
fastq-pipeline/
├── modules/          # Pipeline modules
├── data/            # Your input FASTQ files
├── reference_genome/ # Your reference genome
├── results/         # Pipeline outputs
├── workflow.nf      # Main pipeline script
└── nextflow.config  # Configuration
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

## Contact

GitHub: [@aditijadhav8919-design](https://github.com/aditijadhav8919-design)

## License

MIT License
