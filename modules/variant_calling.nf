/*
 * VARIANT_CALLING - Call variants using bcftools
 */
process VARIANT_CALLING {
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf")
    
    script:
    """
    bcftools mpileup -f ${params.reference} ${bam} | bcftools call -mv -o ${sample_id}.vcf
    """
}
