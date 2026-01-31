// Variant Calling Module

process VARIANT_CALLING {
    
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf
    
    script:
    """
    ${params.bcftools_bin} mpileup -f ${params.ref} ${bam_file} | ${params.bcftools_bin} call -mv -o ${sample_id}.vcf
    """
}
