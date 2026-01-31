#!/usr/bin/env nextflow
/*
 * FASTQ Processing Pipeline
 * Main entry point - imports workflow from workflows/ directory
 */
include { FASTQ_PIPELINE } from './workflows/workflow.nf'

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    FASTQ_PIPELINE(reads_ch)
}
