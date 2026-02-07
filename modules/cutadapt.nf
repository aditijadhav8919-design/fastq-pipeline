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
    ${params.cutadapt_bin} -o trimmed_1.fastq.gz -p trimmed_2.fastq.gz ${reads[0]} ${reads[1]}
}
