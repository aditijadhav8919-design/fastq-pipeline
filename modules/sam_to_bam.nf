/*
 * SAM_TO_BAM - Convert SAM format to compressed BAM format
 */
process SAM_TO_BAM {
    publishDir "${params.outdir}/bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    samtools view -bS ${sam} > ${sample_id}.bam
    """
}
