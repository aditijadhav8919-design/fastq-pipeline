#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC_RAW } from './modules/fastqc_raw.nf'
include { TRIM_ADAPTERS } from './modules/trimming.nf'
include { FASTQC_TRIMMED } from './modules/fastqc_trimmed.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { SAM_TO_BAM } from './modules/sam_to_bam.nf'
include { SORT_BAM } from './modules/sort_bam.nf'
include { VARIANT_CALLING } from './modules/variant_calling.nf'
include { FILTER_VARIANTS } from './modules/filtering.nf'

params.reads = "data/*_R{1,2}.fastq.gz"
params.ref = "${baseDir}/reference_genome/hg38.fa"
params.outdir = "results"

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)
    FASTQC_RAW(reads_ch)
    trimmed_ch = TRIM_ADAPTERS(reads_ch)
    FASTQC_TRIMMED(trimmed_ch)
    aligned_ch = ALIGNMENT(trimmed_ch)
    bam_ch = SAM_TO_BAM(aligned_ch)
    sorted_ch = SORT_BAM(bam_ch)
    variants_ch = VARIANT_CALLING(sorted_ch)
    FILTER_VARIANTS(variants_ch)
}

workflow.onComplete {
    log.info "Pipeline completed successfully!"
}
