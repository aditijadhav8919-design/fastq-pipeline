// Variant Filtering Module
process FILTER_VARIANTS {
    publishDir "${params.outdir}/filtered_variants", mode: "copy"
    
    input:
    tuple val(sample_id), path(vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf"), emit: filtered
    
    script:
    """
    bcftools filter -i 'QUAL>20 && DP>5' ${vcf} > ${sample_id}_filtered.vcf
    """
}
