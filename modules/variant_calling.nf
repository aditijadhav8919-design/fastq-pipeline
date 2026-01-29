// Variant Calling Module
process VARIANT_CALLING {
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf
    
    script:
    """
    bcftools mpileup -f ${params.ref} ${bam} | bcftools call -mv -Ov -o ${sample_id}.vcf
    """
}
