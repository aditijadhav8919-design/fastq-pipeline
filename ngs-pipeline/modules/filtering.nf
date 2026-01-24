process FILTER_VARIANTS {
    tag "Filtering $sample_id"
    publishDir "${params.outdir}/08_filtered", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf)
    
    output:
    path "${sample_id}_filtered.vcf"
    
    script:
    """
    awk '\$6 >= 20 || /^#/' ${vcf} > ${sample_id}_filtered.vcf
    """
}
