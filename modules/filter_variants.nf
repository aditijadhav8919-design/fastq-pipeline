/*
 * FILTER_VARIANTS - Filter variants based on quality
 */
process FILTER_VARIANTS {
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf")
    
    script:
    """
    bcftools filter -i 'QUAL>20' ${vcf} > ${sample_id}.filtered.vcf
    """
}
