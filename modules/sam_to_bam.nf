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
   ${params.samtools_bin} view -bS ${sam_file} > ${sample_id}.bam
    """
}
