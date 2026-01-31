// Mark Duplicates Module

process MARK_DUPLICATES {
    
    publishDir "${params.outdir}/deduplicated", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sorted_bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_dedup.bam"), emit: dedup_bam
    
    script:
    """
    ${params.samtools_bin} markdup -r ${sorted_bam} ${sample_id}_dedup.bam
    """
}
