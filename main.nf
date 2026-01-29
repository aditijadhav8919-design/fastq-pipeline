#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include workflow
include { NGS_PIPELINE } from './workflows/workflow.nf'

// Parameters
params.reads = "test_data/*_R{1,2}.fastq.gz"
params.ref = "${baseDir}/reference_genome/hg38.fa"
params.outdir = "${baseDir}/results"

log.info """
==========================================
 NGS VARIANT CALLING PIPELINE
==========================================
 reads    : ${params.reads}
 reference: ${params.ref}
 outdir   : ${params.outdir}
==========================================
"""

// Main workflow
workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: false)
    NGS_PIPELINE(reads_ch)
}
