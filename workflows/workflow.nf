#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include modules with relative paths
include { FASTQC as FASTQC_RAW } from '../modules/fastqc.nf'
include { CUTADAPT } from '../modules/cutadapt.nf'
include { FASTQC as FASTQC_TRIMMED } from '../modules/fastqc.nf'
include { BWA_ALIGN } from '../modules/bwa_align.nf'
include { SAM_TO_BAM } from '../modules/sam_to_bam.nf'
include { SORT_BAM } from '../modules/sort_bam.nf'
include { VARIANT_CALLING } from '../modules/variant_calling.nf'

workflow QC_PIPELINE {
    
    take:
    reads_ch
    
    main:
    // Quality control on raw reads
    FASTQC_RAW(reads_ch)
    
    // Trim adapters
    CUTADAPT(reads_ch)
    
    // Quality control on trimmed reads
    FASTQC_TRIMMED(CUTADAPT.out)
    
    // Alignment
    BWA_ALIGN(CUTADAPT.out)
    
    // Convert SAM to BAM
    SAM_TO_BAM(BWA_ALIGN.out)
    
    // Sort BAM
    SORT_BAM(SAM_TO_BAM.out)
    
    // Variant calling
    VARIANT_CALLING(SORT_BAM.out)
    
    emit:
    vcf = VARIANT_CALLING.out
}
