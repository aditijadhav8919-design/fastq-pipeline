/*
 * SORT_BAM - Sort BAM files by coordinate
 */
process SORT_BAM {
    publishDir "${params.outdir}/sorted_bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")
    
    script:
    """
    ${params.samtools_bin} sort ${bam_file} -o ${sample_id}_sorted.bam
    """
}
