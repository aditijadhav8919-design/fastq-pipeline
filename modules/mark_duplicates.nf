/*
 * MARK_DUPLICATES - Identify and mark PCR duplicates
 */
process MARK_DUPLICATES {
    publishDir "${params.outdir}/dedup", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam")
    
    script:
    """
    picard MarkDuplicates I=${bam} O=${sample_id}.dedup.bam M=${sample_id}.metrics.txt
    """
}
