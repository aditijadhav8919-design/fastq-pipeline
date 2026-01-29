process SAM_TO_BAM {
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    
    script:
    """
    samtools view -Sb ${sam_file} > ${sample_id}.bam
    """
}
