#!/usr/bin/env nextflow

// Include all modules inside the workflow
include { FASTQC_RAW } from '../modules/fastqc_raw.nf'
include { TRIMMING } from '../modules/trimming.nf'
include { ALIGNMENT } from '../modules/alignment.nf'
include { SORT_BAM } from '../modules/sort_index.nf'
include { VARIANT_CALLING } from '../modules/variant_calling.nf'

workflow FASTQ_PIPELINE {

    // Input channel (reads)
    reads_ch = Channel.fromPath(params.reads)

    // Steps
    qc_raw    = FASTQC_RAW(reads_ch)
    trimmed   = TRIMMING(reads_ch)
    aligned   = ALIGNMENT(trimmed.out)
    sorted    = SORT_BAM(aligned.out)
    variants  = VARIANT_CALLING(sorted.out)
}
