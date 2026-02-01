#!/usr/bin/env nextflow
// Include all modules inside the workflow
include { FASTQC } from '../modules/fastqc.nf'
include { CUTADAPT } from '../modules/cutadapt.nf'
include { BWA_ALIGN } from '../modules/bwa_align.nf'
include { SORT_BAM } from '../modules/sort_bam.nf'
include { VARIANT_CALLING } from '../modules/variant_calling.nf'

workflow FASTQ_PIPELINE {
    take:
    reads_ch
    
    main:
    // Steps
    FASTQC(reads_ch)
    CUTADAPT(reads_ch)
    BWA_ALIGN(CUTADAPT.out)
    SORT_BAM(BWA_ALIGN.out)
    VARIANT_CALLING(SORT_BAM.out)
}
