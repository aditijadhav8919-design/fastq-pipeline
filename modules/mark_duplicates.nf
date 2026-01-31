// Mark Duplicates Module
process MARK_DUPLICATES {
    publishDir "${params.outdir}/deduplicated", mode: "copy"
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}_dedup.bam"), path("${sample_id}_dedup.bam.bai"), emit: dedup
    
    script:
    """
    samtools markdup -r ${bam} ${sample_id}_dedup.bam
    samtools index ${sample_id}_dedup.bam
    """
}
