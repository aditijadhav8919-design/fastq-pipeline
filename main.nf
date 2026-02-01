#!/usr/bin/env nextflow

/*
 * FASTQ Processing Pipeline
 * Main entry point
 */

include { FASTQC_PIPELINE } from './workflows/workflow.nf'

workflow {
    FASTQC_PIPELINE()
}
