#!/usr/bin/env nextflow

// Include all modules with new structure
include { FASTQC } from '../modules/fastqc/main.nf'
include { CUTADAPT } from '../modules/cutadapt/main.nf'
include { SAM_TO_BAM } from '../modules/sam_to_bam/main.nf'
include { SORT_BAM } from '../modules/sort_bam/main.nf'
include { MARK_DUPLICATES } from '../modules/mark_duplicates/main.nf'
include { VARIANT_CALLING } from '../modules/variant_calling/main.nf'
include { FILTER_VARIANTS } from '../modules/filter_variants/main.nf'

workflow FASTQ_PIPELINE {
    take:
    reads_ch

    main:
    // Quality control
    FASTQC(reads_ch)
    
    // Adapter trimming
    CUTADAPT(reads_ch)
    
    // Alignment (you'll need to add a BWA or BOA alignment module)
    // For now, assuming you have aligned SAM files
    
    // Convert SAM to BAM
    SAM_TO_BAM(CUTADAPT.out)
    
    // Sort BAM
    SORT_BAM(SAM_TO_BAM.out)
    
    // Mark duplicates
    MARK_DUPLICATES(SORT_BAM.out)
    
    // Variant calling
    VARIANT_CALLING(MARK_DUPLICATES.out)
    
    // Filter variants
    FILTER_VARIANTS(VARIANT_CALLING.out)
}
