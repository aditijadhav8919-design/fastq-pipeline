process SORT_BAM {
    publishDir "${params.outdir}/sorted", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: sorted
    
    script:
    """
    samtools sort ${bam} -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """
}
