/*
 * CUTADAPT - Trim adapters and low-quality bases from reads
 */
process CUTADAPT {
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz"), emit: trimmed
    
    script:
    
}
