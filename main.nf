#!/usr/bin/env nextflow

/*
 * FASTQ Processing Pipeline
 * Main entry point - imports workflow from workflows/ directory
 */

include { QC_PIPELINE } from './workflows/workflow.nf'

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: false)
    QC_PIPELINE(reads_ch)
}
