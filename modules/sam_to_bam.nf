process SAM_TO_BAM {
    tag "Converting $sample_id"
    publishDir "${params.outdir}/05_bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    samtools view -Sb ${sam} > ${sample_id}.bam
    """
}
