#!/bin/bash

echo "Creating all 8 modules..."

# Module 1: FastQC Raw
cat > modules/fastqc_raw.nf << 'MODULE1'
process FASTQC_RAW {
    tag "QC on $sample_id"
    publishDir "${params.outdir}/01_fastqc_raw", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{html,zip}"
    
    script:
    """
    fastqc -q ${reads}
    """
}
MODULE1

# Module 2: Trimming
cat > modules/trimming.nf << 'MODULE2'
process TRIM_ADAPTERS {
    tag "Trimming $sample_id"
    publishDir "${params.outdir}/02_trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz")
    
    script:
    """
    cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -p ${sample_id}_R2_trimmed.fastq.gz \
        ${reads[0]} ${reads[1]} \
        -m 50 -q 20
    """
}
MODULE2

# Module 3: FastQC Trimmed
cat > modules/fastqc_trimmed.nf << 'MODULE3'
process FASTQC_TRIMMED {
    tag "QC on trimmed $sample_id"
    publishDir "${params.outdir}/03_fastqc_trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*_fastqc.{html,zip}"
    
    script:
    """
    fastqc -q ${reads}
    """
}
MODULE3

# Module 4: Alignment
cat > modules/alignment.nf << 'MODULE4'
process ALIGNMENT {
    tag "Aligning $sample_id"
    publishDir "${params.outdir}/04_aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam")
    
    script:
    """
    bwa mem -t 4 ${params.ref} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}
MODULE4

# Module 5: SAM to BAM
cat > modules/sam_to_bam.nf << 'MODULE5'
process SAM_TO_BAM {
    tag "Converting $sample_id"
    publishDir "${params.outdir}/05_bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    samtools view -Sb ${sam} > ${sample_id}.bam
    """
}
MODULE5

# Module 6: Sort BAM
cat > modules/sort_bam.nf << 'MODULE6'
process SORT_BAM {
    tag "Sorting $sample_id"
    publishDir "${params.outdir}/06_sorted", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")
    
    script:
    """
    samtools sort ${bam} -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """
}
MODULE6

# Module 7: Variant Calling
cat > modules/variant_calling.nf << 'MODULE7'
process VARIANT_CALLING {
    tag "Calling variants $sample_id"
    publishDir "${params.outdir}/07_variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf")
    
    script:
    """
    samtools mpileup -uf ${params.ref} ${bam} > ${sample_id}.bcf
    bcftools view ${sample_id}.bcf > ${sample_id}_raw.vcf
    """
}
MODULE7

# Module 8: Filtering
cat > modules/filtering.nf << 'MODULE8'
process FILTER_VARIANTS {
    tag "Filtering $sample_id"
    publishDir "${params.outdir}/08_filtered", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf)
    
    output:
    path "${sample_id}_filtered.vcf"
    
    script:
    """
    awk '\$6 >= 20 || /^#/' ${vcf} > ${sample_id}_filtered.vcf
    """
}
MODULE8

# Create Main Workflow
cat > workflow.nf << 'WORKFLOW'
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC_RAW } from './modules/fastqc_raw.nf'
include { TRIM_ADAPTERS } from './modules/trimming.nf'
include { FASTQC_TRIMMED } from './modules/fastqc_trimmed.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { SAM_TO_BAM } from './modules/sam_to_bam.nf'
include { SORT_BAM } from './modules/sort_bam.nf'
include { VARIANT_CALLING } from './modules/variant_calling.nf'
include { FILTER_VARIANTS } from './modules/filtering.nf'

params.reads = "data/*_R{1,2}.fastq.gz"
params.ref = "${baseDir}/reference_genome/hg38.fa"
params.outdir = "results"

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)
    FASTQC_RAW(reads_ch)
    trimmed_ch = TRIM_ADAPTERS(reads_ch)
    FASTQC_TRIMMED(trimmed_ch)
    aligned_ch = ALIGNMENT(trimmed_ch)
    bam_ch = SAM_TO_BAM(aligned_ch)
    sorted_ch = SORT_BAM(bam_ch)
    variants_ch = VARIANT_CALLING(sorted_ch)
    FILTER_VARIANTS(variants_ch)
}

workflow.onComplete {
    log.info "Pipeline completed successfully!"
}
WORKFLOW

echo "✅ All modules created!"
echo "📋 Module list:"
ls -lh modules/

echo ""
echo "🔄 Pushing to GitHub..."

# Git commands
git add modules/ workflow.nf
git commit -m "Automated: Added complete 8-step variant calling pipeline"
git push origin main

echo "✅ Successfully pushed to GitHub!"
echo ""
echo "🚀 Ready to run: nextflow run workflow.nf"
