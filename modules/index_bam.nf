process INDEX_BAM {

    input:
    path sorted_bam

    output:
    path "*.bai"

    script:
    """
    samtools index $sorted_bam
    """
}
