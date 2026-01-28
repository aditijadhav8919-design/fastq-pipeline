#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.reads = "data/real_test_R{1,2}.fastq.gz"
params.ref = "${baseDir}/reference_genome/hg38.fa"
params.outdir = "${baseDir}/results"

log.info """
==========================================
 VARIANT CALLING PIPELINE
==========================================
 reads    : ${params.reads}
 reference: ${params.ref}
 outdir   : ${params.outdir}
==========================================
"""

// Process definitions (inline)
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

process FASTQC_TRIMMED {
    tag "QC trimmed $sample_id"
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

// Workflow
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
    log.info """
    ==========================================
    Pipeline ${workflow.success ? 'COMPLETED' : 'FAILED'}
    Duration: ${workflow.duration}
    ==========================================
    """
}
