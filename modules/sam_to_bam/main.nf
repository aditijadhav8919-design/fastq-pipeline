// SAM to BAM conversion Module

process SAM_TO_BAM {
    
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    
    script:
    """
    ${params.samtools_bin} view -bS ${sam_file} > ${sample_id}.bam
    """
}
