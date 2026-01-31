// Variant Filtering Module

process FILTER_VARIANTS {
    
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf"), emit: filtered_vcf
    
    script:
    """
    ${params.bcftools_bin} filter -i 'QUAL>20' ${vcf_file} > ${sample_id}_filtered.vcf
    """
}
