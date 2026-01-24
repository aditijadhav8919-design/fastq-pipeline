#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "data/*_R{1,2}.fastq.gz"
params.ref = "${baseDir}/reference_genome/hg38.fa"
params.outdir = "results"

process fastqc_raw {
    publishDir "${params.outdir}/fastqc_raw"
    input:
    tuple val(sample), path(reads)
    output:
    path "*_fastqc.{zip,html}"
    script:
    "fastqc ${reads}"
}

process align {
    publishDir "${params.outdir}/aligned"
    input:
    tuple val(sample), path(reads)
    output:
    tuple val(sample), path("${sample}.sam")
    script:
    "bwa mem -t 4 ${params.ref} ${reads[0]} ${reads[1]} > ${sample}.sam"
}

process sam_to_bam {
    publishDir "${params.outdir}/bam"
    input:
    tuple val(sample), path(sam)
    output:
    tuple val(sample), path("${sample}.bam")
    script:
    "samtools view -Sb ${sam} > ${sample}.bam"
}

process sort_bam {
    publishDir "${params.outdir}/sorted"
    input:
    tuple val(sample), path(bam)
    output:
    tuple val(sample), path("${sample}_sorted.bam")
    script:
    """
    samtools sort ${bam} -o ${sample}_sorted.bam
    samtools index ${sample}_sorted.bam
    """
}

process call_variants {
    publishDir "${params.outdir}/variants"
    input:
    tuple val(sample), path(bam)
    output:
    tuple val(sample), path("${sample}_variants.txt")
    script:
    """
    samtools mpileup -f ${params.ref} ${bam} > ${sample}_variants.txt
    """
}

process filter_variants {
    publishDir "${params.outdir}/filtered"
    input:
    tuple val(sample), path(variants)
    output:
    path "${sample}_filtered.txt"
    script:
    "awk '\$6 >= 20' ${variants} > ${sample}_filtered.txt"
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)
    
    fastqc_raw(reads_ch)
    aligned = align(reads_ch)
    bam = sam_to_bam(aligned)
    sorted = sort_bam(bam)
    variants = call_variants(sorted)
    filter_variants(variants)
}
