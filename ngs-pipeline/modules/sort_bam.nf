process SORT_BAM {
    tag "Sorting $sample_id"
    publishDir "${params.outdir}/06_sorted", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")
    
    script:
    """
    samtools sort ${bam} -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam
    """
}
