#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC as FASTQC_RAW } from '../modules/fastqc.nf'
include { CUTADAPT } from '../modules/cutadapt.nf'
include { FASTQC as FASTQC_TRIMMED } from '../modules/fastqc.nf'
include { BWA_ALIGN } from '../modules/bwa_align.nf'
include { SAM_TO_BAM } from '../modules/sam_to_bam.nf'
include { SORT_BAM } from '../modules/sort_bam.nf'
include { VARIANT_CALLING } from '../modules/variant_calling.nf'
include { FILTER_VARIANTS } from '../modules/filter_variants.nf'

workflow QC_PIPELINE {
    
    take:
    reads_ch
    
    main:
    FASTQC_RAW(reads_ch)
    CUTADAPT(reads_ch)
    FASTQC_TRIMMED(CUTADAPT.out.trimmed)
    BWA_ALIGN(CUTADAPT.out.trimmed)
    SAM_TO_BAM(BWA_ALIGN.out.sam)
    SORT_BAM(SAM_TO_BAM.out.bam)
    VARIANT_CALLING(SORT_BAM.out.sorted)
    FILTER_VARIANTS(VARIANT_CALLING.out.vcf)
    
    emit:
    vcf = FILTER_VARIANTS.out.filtered
}
