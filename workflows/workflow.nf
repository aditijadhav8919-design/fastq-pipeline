#!/usr/bin/env nextflow

// Include all modules
include { FASTQC } from '../modules/fastqc.nf'
include { CUTADAPT } from '../modules/cutadapt.nf'
include { SAM_TO_BAM } from '../modules/sam_to_bam.nf'
include { SORT_BAM } from '../modules/sort_bam.nf'
include { MARK_DUPLICATES } from '../modules/mark_duplicates.nf'
include { VARIANT_CALLING } from '../modules/variant_calling.nf'
include { FILTER_VARIANTS } from '../modules/filter_variants.nf'
