process VARIANT_CALLING {
    tag "Calling variants $sample_id"
    publishDir "${params.outdir}/07_variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf")
    
    script:
    """
    samtools mpileup -uf ${params.ref} ${bam} > ${sample_id}.bcf
    bcftools view ${sample_id}.bcf > ${sample_id}_raw.vcf
    """
}
