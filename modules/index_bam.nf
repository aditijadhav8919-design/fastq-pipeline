process INDEX_BAM {
    tag "Indexing $sample_id"
    publishDir "${params.outdir}/06_sorted", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sorted_bam)
    
    output:
    tuple val(sample_id), path(sorted_bam), path("${sorted_bam}.bai")
    
    script:
    """
    samtools index ${sorted_bam}
    """
}
